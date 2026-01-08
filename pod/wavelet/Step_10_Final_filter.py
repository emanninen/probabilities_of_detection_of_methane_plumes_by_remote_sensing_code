import json
import cv2
import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import gdal
from shapely.geometry import shape
from scipy.ndimage import binary_dilation
import argparse


def boundary_nan_filter(initial_file, preserved_file, plume_csv_file):
    # Read csv file
    output_df = pd.read_csv(plume_csv_file)
    ratio_nan_within_list = []
    ratio_nan_bound_list = []

    # Read GeoTIFF images
    image_dataset = gdal.Open(initial_file, gdal.GA_ReadOnly)
    image_data = image_dataset.GetRasterBand(1).ReadAsArray()

    preserved_dataset = gdal.Open(preserved_file, gdal.GA_ReadOnly)
    preserved_image = preserved_dataset.GetRasterBand(1).ReadAsArray()
    geo_transform = preserved_dataset.GetGeoTransform()

    # Create connected component labels
    preserved_image_binary = ~np.isnan(preserved_image).astype(bool)
    num_labels, labels = cv2.connectedComponents(preserved_image_binary.astype(np.uint8))

    # Check if the count of clumps in GeoTIFF aligns with the count of rows in CSV
    if num_labels - 1 != len(output_df):
        raise Exception('ERROR: The number of plume masks does not align in GeoTIFF and CSV.')

    for index in range(1, num_labels):
        # Select a clump
        preserved_image_copy = preserved_image.copy()
        preserved_image_copy[labels != index] = np.nan

        # NAN value ratio within the mask
        xch4a_count = np.count_nonzero(~np.isnan(preserved_image_copy))
        preserved_image_copy_2 = np.where(np.isnan(image_data), np.nan, preserved_image_copy)
        xch4a_count_within = np.count_nonzero(~np.isnan(preserved_image_copy_2))

        ratio_nan_within = (xch4a_count - xch4a_count_within) / xch4a_count
        ratio_nan_within_list.append(ratio_nan_within)

        # Create clump mask, expanded clump mask, and expanded boundary mask
        clump_mask = ~np.isnan(preserved_image_copy)

        expansion_size = 1
        structure = np.ones((3, 3))
        expanded_clump_mask = binary_dilation(clump_mask, structure=structure, iterations=expansion_size)
        expanded_area_mask = expanded_clump_mask & ~clump_mask

        # Create expanded boundary images
        bound_preserved_image = np.where(expanded_area_mask, image_data, np.nan)
        bound_label = np.where(expanded_area_mask, labels, np.nan)

        # Calculate NAN boundary pixel count and its count ratio
        count_nan_bound = np.sum(~np.isnan(bound_label)) - np.sum(~np.isnan(bound_preserved_image))
        ratio_nan_bound = count_nan_bound / np.sum(~np.isnan(bound_label))
        ratio_nan_bound_list.append(ratio_nan_bound)

    output_df['Ratio_nan_within'] = ratio_nan_within_list
    output_df['Ratio_nan_bound'] = ratio_nan_bound_list
    output_df.to_csv(plume_csv_file, index=False)

    print('Plume csv file updated: ', plume_csv_file)


def final_filter(plume_csv_file, polygon_json_file, point_json_file, rate_json_file, size_min_threshold, hotspot_size_min_threshold, hotspot_size_ratio_minimum, hotspot_size_ratio_maximum, fibre_length_ratio_threshold, fibre_width_threshold, ratio_nan_within_threshold, ratio_nan_bound_threshold):
    # Define output files
    output_csv_file   = plume_csv_file.split('.csv')[0] +'_v2.csv'
    out_pol_json_file = polygon_json_file.split('.geojson')[0] + '_v2.geojson'
    out_poi_json_file = point_json_file.split('.geojson')[0] + '_v2.geojson'

    # Read input files
    df = pd.read_csv(plume_csv_file)
    print('Number of plumes before filtering: ', len(df))

    pol_gdf = gpd.read_file(polygon_json_file).set_index('Index')
    poi_gdf = gpd.read_file(point_json_file).set_index('plume_id')
    rat_gdf = gpd.read_file(rate_json_file).set_index('plume_id')

    # Merge quantification results to point json by point indices
    poi_gdf['flux'] = rat_gdf.loc[poi_gdf.index, 'flux']
    poi_gdf['flux_sd'] = rat_gdf.loc[poi_gdf.index, 'flux_sd']
    poi_gdf['flux_lo'] = rat_gdf.loc[poi_gdf.index, 'flux_lo']
    poi_gdf['flux_hi'] = rat_gdf.loc[poi_gdf.index, 'flux_hi']

    # Merge quantification results to polygon json by polygon indices
    poi_gdf = poi_gdf.set_index('polygon_index')
    pol_gdf['flux'] = poi_gdf.loc[poi_gdf.index, 'flux']
    pol_gdf['flux_sd'] = poi_gdf.loc[poi_gdf.index, 'flux_sd']
    pol_gdf['flux_lo'] = poi_gdf.loc[poi_gdf.index, 'flux_lo']
    pol_gdf['flux_hi'] = poi_gdf.loc[poi_gdf.index, 'flux_hi']

    # Create new columns in the csv file
    # df['Hotspot_size_ratio'] = df['Size_hotspot'] / df['Size']
    df['Update'] = df['Wind_dir']

    # Filter plumes in the csv file
    df.loc[(df['Fibre_length_ratio'] > fibre_length_ratio_threshold) &
        (df['Fibre_width'] < fibre_width_threshold), 'Update'] = False
    df.loc[df['Hotspot_size_ratio'] < hotspot_size_ratio_minimum, 'Update'] = False
    df.loc[df['Hotspot_size_ratio'] > hotspot_size_ratio_maximum, 'Update'] = True
    df.loc[df['Ratio_nan_within'] > ratio_nan_within_threshold, 'Update'] = False
    df.loc[df['Ratio_nan_bound'] > ratio_nan_bound_threshold, 'Update'] = False
    df.loc[df['Size'] < size_min_threshold, 'Update'] = False
    df.loc[df['Size_hotspot'] < hotspot_size_min_threshold, 'Update'] = False

    # Filter plumes in the geojson file
    pol_gdf = pol_gdf.reset_index()
    poi_gdf = poi_gdf.reset_index()

    pol_gdf[['Wind_dir', 'Update']] = df[['Wind_dir', 'Update']]
    pol_gdf['Pass_NAN_filter'] = df['Ratio_nan_bound'] <= ratio_nan_bound_threshold

    filtered_df = df[df['Update'] == True]
    filtered_pol_gdf = pol_gdf[pol_gdf['Index'].isin(filtered_df['Index'])]
    filtered_poi_gdf = poi_gdf[poi_gdf['polygon_index'].isin(filtered_df['Index'])]
    print('Number of plumes after filtering:  ', len(filtered_df))

    # Save output
    df.to_csv(output_csv_file, index=False)
    filtered_pol_gdf.to_file(out_pol_json_file, driver='GeoJSON')
    filtered_poi_gdf.to_file(out_poi_json_file, driver='GeoJSON')

    print('Plume csv file updated:       ', output_csv_file)
    print('Polygon GeoJSON file created: ', out_pol_json_file)
    print('Point GeoJSON file created:   ', out_poi_json_file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Final filter.")
    parser.add_argument("flight", type=str, help="Flight name.")
    
    args = parser.parse_args()
    flight = args.flight

    initial_file =      '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' +  flight + '.tif'
    preserved_file =    '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' +  flight + '_preserved_v2.tif'
    plume_csv_file =    '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered_plumes.csv'
    polygon_json_file = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered.geojson'
    point_json_file   = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered_plumes.geojson'
    rate_json_file    = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + 'growing_boxes.geojson'
    
    # Determine filtering threshold values
    size_min_threshold =           200
    hotspot_size_min_threshold =   15
    hotspot_size_ratio_minimum =   0.003
    hotspot_size_ratio_maximum =   0.03
    fibre_length_ratio_threshold = 1.25
    fibre_width_threshold =        0.0025
    ratio_nan_within_threshold =   0.5
    ratio_nan_bound_threshold =    0.25
    
    
    # Step 1: Run boundary nan filter
    boundary_nan_filter(initial_file, preserved_file, plume_csv_file)

    # Step 2: Run final filter
    final_filter(plume_csv_file, polygon_json_file, point_json_file, rate_json_file, 
    size_min_threshold, hotspot_size_min_threshold, hotspot_size_ratio_minimum, 
    hotspot_size_ratio_maximum, fibre_length_ratio_threshold, fibre_width_threshold, 
    ratio_nan_within_threshold, ratio_nan_bound_threshold)
    