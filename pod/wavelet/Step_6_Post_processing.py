import os
import numpy as np
import pandas as pd
import numpy.ma as ma
from osgeo import gdal, gdal_array, ogr, osr
import cv2
import geojson
import json
from shapely.geometry import shape, Point, Polygon
import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_mask
from skimage.morphology import skeletonize
import networkx
from shapely.geometry import LineString
from Image_proc_functions import *
from Detached_plumes_copy import *
import argparse


def post_processing(initial_file, preserved_file, point_jsonfile, base_hotspot_percentile, threshold_distance = 0.001):
    # Read GeoTIFF image
    image_dataset = gdal.Open(initial_file, gdal.GA_ReadOnly)
    image_data = image_dataset.GetRasterBand(1).ReadAsArray()

    preserved_dataset = gdal.Open(preserved_file, gdal.GA_ReadOnly)
    preserved_image = cv2.imread(preserved_file, cv2.IMREAD_UNCHANGED)

    # Get geotransform info
    cols = preserved_dataset.RasterXSize
    rows = preserved_dataset.RasterYSize
    geo_transform = preserved_dataset.GetGeoTransform()
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
    dat2 = np.stack((lon_array, lat_array, preserved_image), axis=-1)

    preserved_image_binary = ~np.isnan(preserved_image).astype(bool)
    num_labels, labels = cv2.connectedComponents(preserved_image_binary.astype(np.uint8))

    # Calculate base_hotspot_threshold
    base_hotspot_threshold = np.nanpercentile(image_data, base_hotspot_percentile)
    print('Base_hotspot_threshold:            ', base_hotspot_threshold)

    base_background_mean = np.nanmean(image_data)
    print('Mean background mean:              ', base_background_mean)

    base_background_std = np.nanstd(image_data)
    print('Mean background std:               ', base_background_std)
    
    # Read point GeoJSON data
    with open(point_jsonfile, 'r') as f:
        points_data = json.load(f)

    flight_name = point_jsonfile.split('/')[-3]
    
    # Create an output masked image
    output_masked_image = np.zeros_like(preserved_image)
    output_masked_image[output_masked_image == 0] = np.nan
    
    # Create polygons from clumps in the masked image
    preserved_image_binary = ~np.isnan(preserved_image).astype(bool)
    num_labels, labels = cv2.connectedComponents(preserved_image_binary.astype(np.uint8))
    
    polygons_data = labeled_image_to_geojson(num_labels, labels, geo_transform)

    # Extract polygons and points
    polygons = [shape(feature['geometry']) for feature in polygons_data['features']]
    points = [shape(feature['geometry']) for feature in points_data['features']]
    
    # Loop through all the polygons
    filtered_polygons = []
    results = []
    new_index = 0
    for index, polygon in enumerate(polygons):
        # Find the closest point to this polygon as its point source
        min_dist, closest_point = closest_point_to_polygon(polygon, points)
        if min_dist <= threshold_distance:
            # The distance between polygon and point should be < threshold_distance
            filtered_polygons.append(polygon)
            
            # Preserve this single mask
            preserved_image_copy = preserved_image.copy()
            preserved_image_copy[labels != index + 1] = np.nan

            # Add this mask to the output masked image
            output_masked_image = np.where(~np.isnan(preserved_image_copy), preserved_image_copy, output_masked_image)
            new_index += 1

            # Size (count of non-nan values)
            xch4a_count = np.count_nonzero(~np.isnan(preserved_image_copy))
            
            # Find its centroid
            center_lon = np.mean(np.ma.masked_array(dat2[:, :, 0], mask=np.isnan(preserved_image_copy)))
            center_lat = np.mean(np.ma.masked_array(dat2[:, :, 1], mask=np.isnan(preserved_image_copy)))

            # Crop the image data to a sub image
            half_box_size = 0.01
            sub_original_image, start_x, start_y, width, height = crop_geotiff_image(image_data, 
                                                                                    initial_file, center_lon, center_lat, half_box_size)

            # Preserve only the sub image background surrounding the mask
            x_s, x_e, y_s, y_e = start_y, start_y+height, start_x, start_x+width
            sub_labels = labels[x_s:x_e, y_s:y_e]
            sub_original_image_copy = np.where(sub_labels == 0, sub_original_image, np.nan)
            sub_rows, sub_cols = sub_original_image_copy.shape
            
            # Check background pixel count, enlarge the sub image if background pixels too few
            background_count = np.count_nonzero(~np.isnan(sub_original_image_copy))
            background_count_ratio = background_count / (sub_rows * sub_cols)
            
            # Enlarge the sub image if background pixel count < half of the total pixel count
            # and background pixel count < large_label_mask_threshold
            large_label_mask_threshold = 30000

            while background_count_ratio < 0.5 and background_count < large_label_mask_threshold:
                half_box_size += 0.005
                sub_original_image, start_x, start_y, width, height = crop_geotiff_image(image_data, 
                                                                                    initial_file, center_lon, center_lat, half_box_size)

                # Preserve only the sub image background surrounding the mask
                x_s, x_e, y_s, y_e = start_y, start_y+height, start_x, start_x+width
                sub_labels = labels[x_s:x_e, y_s:y_e]
                sub_original_image_copy = np.where(sub_labels == 0, sub_original_image, np.nan)
                sub_rows, sub_cols = sub_original_image_copy.shape
                
                # Check background pixel count
                background_count = np.count_nonzero(~np.isnan(sub_original_image_copy))
                background_count_ratio = background_count / (sub_rows * sub_cols)
            
            # Calculate background std
            sub_background_std = np.nanstd(sub_original_image_copy)

            # Hotspot threshold and hotspot size
            # Gaussian distribution: mu + 3.4 * sigma = 99.97th percentile
            std_diff_threshold = 2
            std_to_percentile_scale = 3.4

            if (sub_background_std - base_background_std) > std_diff_threshold:
                hotspot_threshold = base_hotspot_threshold + std_to_percentile_scale * (sub_background_std - base_background_std)
                hotspot_count = np.nansum(preserved_image_copy > hotspot_threshold)
            else:
                hotspot_threshold = base_hotspot_threshold
                hotspot_count = np.nansum(preserved_image_copy > base_hotspot_threshold)
            
            # Hotspot size ratio
            hotspot_size_ratio = hotspot_count / xch4a_count

            # Fibre length ratio and fibre width
            fibre_length, _ = calculate_fibre_length(polygon)
            major_axis_length, _ = get_major_axis(polygon)

            fibre_to_major_ratio = fibre_length / major_axis_length
            polygon_area = polygon.area
            fibre_width = polygon_area / fibre_length

            # Export key info to a DataFrame
            results.append({'Index': new_index, 'Latitude': closest_point.y, 'Longitude': closest_point.x, 'Hotspot_thrld': hotspot_threshold,
                            'Sub_bck_std': sub_background_std, 'Size': xch4a_count, 'Size_hotspot': hotspot_count, 
                            'Hotspot_size_ratio': hotspot_size_ratio, 'Fibre_length_ratio':fibre_to_major_ratio, 'Fibre_width': fibre_width})
            
    # Save outputs
    output_masked_file = preserved_file.split('.tif')[0] + '_v2.tif'
    csv_file = point_jsonfile.split('/plumes.geojson')[0] +  '/' + flight_name + '_filtered_plumes.csv'
    filtered_jsonfile = point_jsonfile.split('/plumes.geojson')[0] +  '/' + flight_name + '_filtered.geojson'
    
    # Save output masked image
    create_geotiff(preserved_dataset, output_masked_image, output_masked_file)

    # Write results to a CSV file
    df = pd.DataFrame(results)
    df.to_csv(csv_file, index=False)
    print('Plume csv file created:       ', csv_file)

    # Write filtered polygons to a GeoJSON file
    filtered_features = [{'type': 'Feature', 'geometry': polygon.__geo_interface__, 'properties': {
                'Index': result['Index'],
                'Latitude': result['Latitude'],
                'Longitude': result['Longitude'],
                'Hotspot_thrld': result['Hotspot_thrld'],
                'Sub_bck_std': float(result['Sub_bck_std']),
                'Size': int(result['Size']),
                'Size_hotspot': int(result['Size_hotspot']),
                'Hotspot_size_ratio': result['Hotspot_size_ratio'],
                'Fibre_length_ratio': result['Fibre_length_ratio'],
                'Fibre_width': result['Fibre_width'],
            }
        } for polygon, result in zip(filtered_polygons, results)]
    filtered_geojson = {'type': 'FeatureCollection', 'features': filtered_features}

    with open(filtered_jsonfile, 'w') as f:
        json.dump(filtered_geojson, f)
    print('Polygon GeoJSON file created: ', filtered_jsonfile)


def calculate_fibre_length(polygon, scale=10000):  
    # Convert polygon to binary image
    minx, miny, maxx, maxy = polygon.bounds
    width = int((maxx - minx) * scale)
    height = int((maxy - miny) * scale)
    
    img = np.zeros((height, width), dtype=np.uint8)
    coords = np.array([(int((x - minx) * scale), int((maxy - y) * scale)) for x, y in polygon.exterior.coords], dtype=np.int32)
    cv2.fillPoly(img, [coords], 1)

    # Skeletonize the binary image
    skeleton = skeletonize(img).astype(np.uint8)

    # Extract coordinates of skeleton pixels
    y_coords, x_coords = np.nonzero(skeleton)
    coords = list(zip(x_coords, y_coords))

    # Create a graph from the skeleton pixels
    G = networkx.Graph()
    for x, y in coords:
        G.add_node((x, y))
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                nx, ny = x + dx, y + dy
                if (nx, ny) in G.nodes:
                    G.add_edge((x, y), (nx, ny))

    # Find the longest path in the skeleton
    lengths = dict(networkx.all_pairs_shortest_path_length(G))
    max_length = 0
    max_path = None
    for node in lengths:
        for target in lengths[node]:
            if lengths[node][target] > max_length:
                max_length = lengths[node][target]
                max_path = networkx.shortest_path(G, source=node, target=target)

    if max_path is None:
        return None, None

    # Convert path to coordinates in the original scale
    fibre_line = [(x / scale + minx, maxy - y / scale) for x, y in max_path]
    
    # Scale the length back to the original size
    fibre_length = max_length / scale

    return fibre_length, fibre_line


def get_major_axis(polygon):
    # minx, miny, maxx, maxy = polygon.bounds
    longest_line = LineString([polygon.exterior.coords[0], polygon.exterior.coords[2]])
    
    for i in range(len(polygon.exterior.coords)):
        for j in range(i + 1, len(polygon.exterior.coords)):
            line = LineString([polygon.exterior.coords[i], polygon.exterior.coords[j]])
            if line.length > longest_line.length:
                longest_line = line
    
    return longest_line.length, longest_line


def load_polygons_from_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    polygons = [shape(feature['geometry']) for feature in data['features'] if feature['geometry']['type'] == 'Polygon']
    return polygons





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DI-derived plume origin distance filtering.")
    parser.add_argument("flight", type=str, help="Flight name.")
    parser.add_argument("base_hotspot_percentile", nargs='?', type=float, default = 99.97, 
                        help="Base hotspot percentile threhold")
    
    # Run the function with the provided arguments
    args = parser.parse_args()
    flight = args.flight
    base_hotspot_percentile = args.base_hotspot_percentile


    initial_file =       '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '.tif'
    preserved_file     = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '_preserved.tif'
    point_jsonfile     = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/plumes.geojson'
    # initial_file =       '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '.tif'
    # preserved_file =     '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_preserved_v2_v7.tif'
    # point_jsonfile =     '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/plumes_test.geojson'


    post_processing(initial_file, preserved_file, point_jsonfile, base_hotspot_percentile)
