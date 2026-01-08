#%%
import os
import os.path
import math
import numpy as np
import numpy.ma as ma
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_mask
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from Image_proc_functions import *
from Wavelet_functions import *
import cv2
import json
from shapely.geometry import Point
import argparse
import warnings
warnings.filterwarnings("ignore")


def mask_processing(initial_file, masked_file, large_label_mask_threshold, scale, sigma_factor):
    # Read GeoTIFF files
    image_dataset = gdal.Open(initial_file, gdal.GA_ReadOnly)
    image_data = image_dataset.GetRasterBand(1).ReadAsArray()
    masked_dataset = gdal.Open(masked_file, gdal.GA_ReadOnly)
    masked_image = masked_dataset.GetRasterBand(1).ReadAsArray()

    # Get geotransform info
    cols = masked_dataset.RasterXSize
    rows = masked_dataset.RasterYSize
    geo_transform = masked_dataset.GetGeoTransform()
    min_x = geo_transform[0]
    max_y = geo_transform[3]

    # Create variable dat2 [lon, lat, xch4a]
    lon_array = np.zeros((rows, cols))
    lat_array = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            lon, lat = gdal.ApplyGeoTransform(geo_transform, j, i)
            lon_array[i, j] = lon
            lat_array[i, j] = lat
    dat2 = np.stack((lon_array, lat_array, masked_image), axis=-1)

    # Create a hotspot_coeff image
    hotspot_coeff_image = np.zeros_like(masked_image)

    # Create an output masked image
    output_masked_image = masked_image.copy()

    # Check the size of the largest mask
    masked_image_binary = ~np.isnan(output_masked_image).astype(bool)
    _, labels_to_check = cv2.connectedComponents(masked_image_binary.astype(np.uint8))
    sizes = np.sort(np.bincount(labels_to_check.ravel()))

    while sizes[-2] > large_label_mask_threshold:
        # Preserve the large masks only
        large_masked_image = remove_small_clumps(output_masked_image, large_label_mask_threshold)
        masked_image_binary = ~np.isnan(large_masked_image).astype(bool)
        num_labels, labels = cv2.connectedComponents(masked_image_binary.astype(np.uint8))

        # plt.figure(figsize = (7,7))
        # plt.imshow(labels)
        # plt.colorbar()
        # plt.show()
        
        # Create polygons for the masks
        polygons_data = labeled_image_to_geojson(num_labels, labels, geo_transform)
        polygons = [shape(feature['geometry']) for feature in polygons_data['features']]
        
        # Update hotspot_coeff_image
        hotspot_coeff_image += masked_image_binary.astype(np.uint8)

        # plt.figure(figsize = (7,7))
        # plt.imshow(hotspot_coeff_image)
        # plt.colorbar()
        # plt.show()

        for i in range(num_labels - 1):
            # Preserve a single large mask
            large_masked_image_copy = large_masked_image.copy()
            large_masked_image_copy[labels != i + 1] = np.nan

            # Find its centroid
            center_lon = np.mean(np.ma.masked_array(dat2[:, :, 0], mask=np.isnan(large_masked_image_copy)))
            center_lat = np.mean(np.ma.masked_array(dat2[:, :, 1], mask=np.isnan(large_masked_image_copy)))
            centroid_point = Point(center_lon, center_lat)

            # Use centroid and polygon vertice to find cropping box size
            polygon = polygons[i]
            half_box_size = farthest_dist_point_polygon(polygon, centroid_point) + 0.001

            # Crop the image data to a sub image
            sub_original_image, start_x, start_y, width, height = crop_geotiff_image(
            image_data, initial_file, center_lon, center_lat, half_box_size)
            
            # Run regional wavelet to the cropped original image to get new masks
            min_x_crop = center_lon - half_box_size
            max_y_crop = center_lat + half_box_size
            new_sub_masked_image = regional_wavelet(sub_original_image, initial_file, min_x_crop, max_y_crop, scale, sigma_factor)

            # Replace sub original masks with new masks
            temp_image = output_masked_image.copy()
            output_masked_image[start_y:start_y+height, start_x:start_x+width] = new_sub_masked_image
            output_masked_image[labels != i + 1] = temp_image[labels != i + 1]
        
        # plt.figure(figsize = (7,7))
        # plt.imshow(output_masked_image)
        # plt.colorbar()
        # plt.show()
        
        # Update sizes
        masked_image_binary = ~np.isnan(output_masked_image).astype(bool)
        _, labels_to_check = cv2.connectedComponents(masked_image_binary.astype(np.uint8))
        sizes = np.sort(np.bincount(labels_to_check.ravel()))
    
    # Remove small clumps and fill NAN holes within masks
    label_mask_threshold = 100
    output_masked_image = remove_small_clumps(output_masked_image, label_mask_threshold)

    masked_image_geotransform = (min_x, max_y, rows, cols)
    output_masked_image = fill_nan_within_mask(output_masked_image, masked_image_geotransform, initial_file)
    
    # Update hotspot_coeff_image
    masked_image_binary = ~np.isnan(output_masked_image).astype(bool)
    hotspot_coeff_image += masked_image_binary.astype(np.uint8)

    # Save output masked image
    update_masked_file = masked_file.split('.tif')[0] + '_v2.tif'
    create_geotiff(masked_dataset, output_masked_image, update_masked_file)

    # hotspot_coeff_file = masked_file.split('.tif')[0] + '_hotspot_coeff.tif'
    # hotspot_coeff_image[np.isnan(image_data)] = np.nan
    # create_geotiff(image_dataset, hotspot_coeff_image, hotspot_coeff_file)


def regional_wavelet(image_data, initial_file, min_x, max_y, scale, sigma_factor):
    # Fill NANs in the original image data
    non_nan_values = image_data[~np.isnan(image_data)]
    hist, bins = np.histogram(non_nan_values, bins=50)
    peak_value = bins[np.argmax(hist)]
    noise = np.random.normal(0, np.std(non_nan_values), image_data.shape)
    image_data = np.where(np.isnan(image_data), peak_value + noise, image_data)
    rows, cols = image_data.shape

    # Define wavelet parameters
    wavelet = 'haar'
    denoise_wavelet = 'bior3.5'
    n = math.floor(np.log2(rows)) // 2
    noiseSigma = 16  # factor of the denoising soft threshold (range: 16, 32, 64, a higher value helps remove more small masks that are likely to be noise)
    label_mask_threshold = 100  # = 10p*10p = 100m*100m (avoid high values: may miss small-shape plumes)
    gaussian_sigma = 3  # sigma of Gaussian filtering (avoid low values (<1): masks may be too "shattered")

    # Image processing using wavelet transform
    _, _, _, denoised_image = wavelet_processing(image_data, wavelet, 
                                                 denoise_wavelet, scale, n, noiseSigma, (1, 0, 0, 0))

    # Add Gaussian filter to denoised image
    denoised_image = gaussian_filter(denoised_image, sigma=gaussian_sigma)

    # Create masks of original image based on denoised image
    masked_image = plume_masking(denoised_image, image_data, sigma_factor, label_mask_threshold)
    masked_image = masked_image.filled(np.nan)
    masked_image[np.isnan(image_data)] = np.nan
    
    # Further mask out clumps that are smaller than threshold
    masked_image = remove_small_clumps(masked_image, label_mask_threshold)

    # plt.figure(figsize = (7,7))
    # plt.imshow(masked_image)
    # plt.colorbar()
    # plt.show()

    # Fill NAN holes within the masks
    if np.any(~np.isnan(masked_image)):
        masked_image_geotransform = (min_x, max_y, rows, cols)
        filled_nan_masked_image = fill_nan_within_mask(masked_image, masked_image_geotransform, initial_file)
    else:
        filled_nan_masked_image = masked_image

    return filled_nan_masked_image





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mask processing (regional wavelet).")
    parser.add_argument("flight", type=str, help="Flight name.")
    parser.add_argument("scale", type=float, help="Factor of the consistent signal threshold")
    parser.add_argument("sigma_factor", type=float, help="Factor of the plume masking threshold")
    parser.add_argument("large_label_mask_threshold", nargs='?', type=int, default = 30000, 
                        help="Maximum mask size threshold")
    
    # Run the function with the provided arguments
    args = parser.parse_args()
    flight = args.flight
    scale = args.scale
    sigma_factor = args.sigma_factor
    large_label_mask_threshold = args.large_label_mask_threshold
    
    initial_file = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '.tif'
    masked_file =  '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '_masked.tif'
    # initial_file = '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '.tif'
    # masked_file =  '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_masked_test_1.5.tif'
    

    mask_processing(initial_file, masked_file, large_label_mask_threshold, scale, sigma_factor)




# flight = 'MSAT_240813'
# scale, sigma_factor = 2, 1.5
# initial_file = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '.tif'
# masked_file =  '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '_masked.tif'
# large_label_mask_threshold = 6666
# mask_processing(initial_file, masked_file, large_label_mask_threshold, scale, sigma_factor)
# %%
