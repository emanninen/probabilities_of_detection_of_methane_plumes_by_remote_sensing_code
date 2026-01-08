import os
import numpy as np
from osgeo import gdal, gdal_array, ogr, osr
import cv2
from shapely.geometry import shape, Point, Polygon
import matplotlib.pyplot as plt
from Image_proc_functions import *
import argparse


def remove_mask_base_hotspot_threshold(
    initial_file
    , masked_file
    , base_hotspot_percentile
    , hotspot_label_mask_threshold
    , base_hotspot_threshold = 0
):
    # Read GeoTIFF images
    image_dataset = gdal.Open(initial_file, gdal.GA_ReadOnly)
    image_data = image_dataset.GetRasterBand(1).ReadAsArray()

    masked_dataset = gdal.Open(masked_file, gdal.GA_ReadOnly)
    masked_image = masked_dataset.GetRasterBand(1).ReadAsArray()

    if base_hotspot_threshold == 0:
        # Calculate base_hotspot_threshold
        base_hotspot_threshold = np.nanpercentile(image_data, base_hotspot_percentile)

    print('Base_hotspot_threshold:            ', base_hotspot_threshold)


    # Remove masks based on base_hotspot_threshold
    hotspot_image = hotpost_masking(masked_image, base_hotspot_threshold, hotspot_label_mask_threshold)
    preserved_image = preserve_hotspot(masked_image, hotspot_image)
    
    # Save output image
    output_mask_file = masked_file.split('_v2.tif')[0] + '_v3.tif'
    create_geotiff(masked_dataset, preserved_image, output_mask_file)


# Input masked_file: the updated masked file after remove_mask_base_hotspot_threshold
def remove_mask_hotspot_threshold(
    initial_file
    , masked_file
    , base_hotspot_percentile
    , hotspot_label_mask_threshold
    , base_hotspot_threshold = 0
):
    # Read GeoTIFF images
    image_dataset = gdal.Open(initial_file, gdal.GA_ReadOnly)
    image_data = image_dataset.GetRasterBand(1).ReadAsArray()

    masked_dataset = gdal.Open(masked_file, gdal.GA_ReadOnly)
    masked_image = masked_dataset.GetRasterBand(1).ReadAsArray()

    # Get geotransform info
    cols = masked_dataset.RasterXSize
    rows = masked_dataset.RasterYSize
    geo_transform = masked_dataset.GetGeoTransform()

    # Create variable dat2 [lon, lat, xch4a]
    lon_array = np.zeros((rows, cols))
    lat_array = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            lon, lat = gdal.ApplyGeoTransform(geo_transform, j, i)
            lon_array[i, j] = lon
            lat_array[i, j] = lat
    dat2 = np.stack((lon_array, lat_array, masked_image), axis=-1)

    # Calculate base thresholds
    if base_hotspot_threshold ==0:

        base_hotspot_threshold = np.nanpercentile(image_data, base_hotspot_percentile)

    print('Base_hotspot_threshold:            ', base_hotspot_threshold)

    base_background_std = np.nanstd(image_data)
    
    # Define an output hotspot image
    output_hotspot_image = np.zeros_like(masked_image)
    output_hotspot_image[output_hotspot_image == 0] = np.nan
    
    # Create mask clump labels
    masked_image_binary = ~np.isnan(masked_image).astype(bool)
    num_labels, labels = cv2.connectedComponents(masked_image_binary.astype(np.uint8))
    
    # Loop through all the masks
    for i in range(num_labels - 1):
        test_label = i + 1
            
        # Preserve a single mask
        masked_image_copy = masked_image.copy()
        masked_image_copy[labels != test_label] = np.nan

        ### base hotspot thresh
        hotspot_threshold = base_hotspot_threshold

        # Create hotspot image of this single mask
        hotspot_image = hotpost_masking(masked_image_copy, hotspot_threshold, hotspot_label_mask_threshold)
        
        # Add the hotspot of this single mask to the output hotspot image
        output_hotspot_image = np.where(~np.isnan(hotspot_image), hotspot_image, output_hotspot_image)
    
    output_hotspot_image[output_hotspot_image == 0] = np.nan
    output_hotspot_image = remove_small_clumps(output_hotspot_image, hotspot_label_mask_threshold)

    return output_hotspot_image


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hotspot filtering.")
    parser.add_argument("flight", type=str, help="Flight name.")
    parser.add_argument("label_mask_threshold", type=int, 
                        help="Minimum mask size threshold, usually 100 or 200")
    parser.add_argument("hotspot_label_mask_threshold", type=int, 
                        help="Minimum hotspot size threshold, usually 10 or 15")
    parser.add_argument("base_hotspot_percentile", nargs='?', type=float, default = 99.97, help="Base hotspot percentile threhold")
    parser.add_argument("base_hotspot_threshold", nargs='?', type=float, default = 0, help="Base hotspot threhold")
    parser.add_argument("initial_file", nargs='?', type=str, help="input to Step 1")
    parser.add_argument("output_masked_file", nargs='?', type=str, 
                        help="out put of Step 1")
    parser.add_argument("output_dir", nargs='?', type=str, 
                        help="output dir")
    
    # Run the function with the provided arguments
    args = parser.parse_args()
    initial_file = args.initial_file
    output_masked_file = args.output_masked_file
    output_dir = args.output_dir
    flight = args.flight
    label_mask_threshold = args.label_mask_threshold
    hotspot_label_mask_threshold = args.hotspot_label_mask_threshold
    base_hotspot_percentile = args.base_hotspot_percentile
    base_hotspot_threshold = args.base_hotspot_threshold
    
    # output_masked_file = output_dir + flight + '_masked_v3.tif'
    preserved_file =     output_dir + flight + '_preserved.tif'
    # initial_file =       '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '.tif'
    # masked_file =        '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_masked_v5.tif'
    # output_masked_file = '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_masked_v6.tif'
    # preserved_file =     '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_preserved_v2_v7.tif'
    
    # Print base threholds
    image_dataset = gdal.Open(initial_file, gdal.GA_ReadOnly)
    image_data = image_dataset.GetRasterBand(1).ReadAsArray()

    # base_hotspot_threshold = np.nanpercentile(image_data, base_hotspot_percentile)
    print('Base hotspot threshold:            ', base_hotspot_threshold)

    base_background_mean = np.nanmean(image_data)
    print('Mean background mean:              ', base_background_mean)

    base_background_std = np.nanstd(image_data)
    print('Mean background std:               ', base_background_std)
    
    # # Step 1: Remove masks by base hotspot threholds
    # remove_mask_base_hotspot_threshold(initial_file, masked_file, base_hotspot_percentile, hotspot_label_mask_threshold)

    # initial_file
    # , masked_file
    # , base_hotspot_percentile
    # , hotspot_label_mask_threshold
    # , base_hotspot_threshold = 0
    # Step 2: Remove masks by updated hotspot threholds
    output_hotspot_image = remove_mask_hotspot_threshold(
        initial_file = initial_file
        , masked_file = output_masked_file
        , base_hotspot_percentile = base_hotspot_percentile
        , hotspot_label_mask_threshold = hotspot_label_mask_threshold
        , base_hotspot_threshold = base_hotspot_threshold
    )#end remove mask hotspot

    # Step 3: Create a new masked image with preserved masks
    masked_dataset = gdal.Open(output_masked_file, gdal.GA_ReadOnly)
    masked_image = masked_dataset.GetRasterBand(1).ReadAsArray()

    preserved_image = preserve_hotspot(masked_image, output_hotspot_image)
    preserved_image = remove_small_clumps(preserved_image, label_mask_threshold)
    create_geotiff(masked_dataset, preserved_image, preserved_file)



