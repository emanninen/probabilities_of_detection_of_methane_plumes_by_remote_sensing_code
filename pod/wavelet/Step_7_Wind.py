import os
import numpy as np
import pandas as pd
import geopandas as gpd
from netCDF4 import Dataset
from msat_level4_divergence_integral.io.hrrr_input import extract_met_data
from msat_level4_divergence_integral.io.layer_averaging import select_layer
from msat_level4_divergence_integral.divergence_integral.integration import prep_met, _prep_ch4
from msat_level4_divergence_integral.io.l3_input import format_level3
from msat_level4_divergence_integral.plume_finder.plume_characterization import _wind_angle
from msat_level4_divergence_integral.divergence_integral.types import DivergenceIntegralOptions
import json
from shapely.geometry import shape, Polygon
from skimage.draw import polygon as draw_polygon
from skimage.measure import regionprops, label
import matplotlib.pyplot as plt
import argparse

def create_plume_wind_data(directory_path, polygon_file):
    # Find centroid coord for L3 clumps
    gdf = gpd.read_file(polygon_file)
    centroids = gdf['geometry'].centroid
    lonlat_centroids = np.array([(centroid.x, centroid.y) for centroid in centroids])
    print('Number of clumps in the polygon geojson: ', len(lonlat_centroids))

    # Create an output DataFrame
    indices = gdf['Index'].to_numpy()
    df = pd.DataFrame(indices, columns = ['Index'])

    # Loop through all .grib2 files
    files = os.listdir(directory_path)
    for file in files:
        if file.endswith('.grib2'):
            hour = file.split('z.wrfnatf.grib2')[0][-2:]
            hrrr_path = os.path.join(directory_path, file)

            # Load inputs: HRRR data
            met_data = extract_met_data(hrrr_path)
            layer_met_data = select_layer(met_data)
            prepared_met_data = prep_met(layer_met_data)

            wind_interpolator = prepared_met_data.wind_interpolator

            # # Match HRRR with L3 to get wind data X*Y*2[u,v]
            # wind = wind_interpolator(ch4_vertices)

            # Calculate wind angles for L3 clumps
            met_wind_angles = []
            for centroid_lonlat in lonlat_centroids:
                met_wind_angle = _wind_angle(centroid_lonlat, wind_interpolator)[0]
                met_wind_angles.append(met_wind_angle)

            df[hour] = met_wind_angles
    
    output_wind_csv_file = directory_path + '/' + polygon_file.split('/')[-2] + '/wind.csv'
    df.to_csv(output_wind_csv_file, index=False)
    
    print('Wind csv file created: ', output_wind_csv_file)






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract HRRR wind data.")
    parser.add_argument("flight", type=str, help="Flight name.")
    
    args = parser.parse_args()
    flight = args.flight

    directory_path =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/'
    polygon_file =     '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered.geojson'
    wind_file =        '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' +'wind.csv'
    # directory_path =   '/n/home12/zhanzhang/Test_DI/' + flight + '/'
    # polygon_file =     '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' + flight + '_filtered_test.geojson'
    # wind_file =        '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' +'wind.csv'

    create_plume_wind_data(directory_path, polygon_file)
