#%%
import os
import os.path
import math
import json
import numpy as np
import numpy.ma as ma
import pandas as pd
import cv2
from osgeo import gdal, gdal_array, ogr, osr
from shapely.geometry import shape, Point, Polygon, mapping
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from Detached_plumes_copy import *
from Image_proc_functions import *
import argparse

def wind_angle_align(wind_array, angle, degree_incre):
    direction = False
    if (np.max(wind_array) - np.min(wind_array)) <= 180:
        wind_range_min = np.min(wind_array) - degree_incre
        wind_range_max = np.max(wind_array) + degree_incre

        if (wind_range_min >= 0 and wind_range_max <= 360) and (wind_range_min <= angle <= wind_range_max):
            print('Pass wind direction filter.')
            print('HRRR wind angle range: ', (wind_range_min, wind_range_max))
            print('plume wind angle:      ', angle)
            direction = True

        elif (wind_range_min < 0 and wind_range_max <= 360) and ((0 <= angle <= wind_range_max) 
                                                                 or ((wind_range_min % 360) <= angle <= 360)):
            print('Pass wind direction filter.')
            print('HRRR wind angle range: ', (wind_range_min, wind_range_max))
            print('plume wind angle:      ', angle)
            direction = True

        elif (wind_range_min >= 0 and wind_range_max > 360) and ((wind_range_min <= angle <= 360) 
                                                                 or (0 <= angle <= (wind_range_max % 360))):
            print('Pass wind direction filter.')
            print('HRRR wind angle range: ', (wind_range_min, wind_range_max))
            print('plume wind angle:      ', angle)
            direction = True
        else:
            print('Fail wind direction filter')
            wind_range_min = (wind_range_min % 360) if wind_range_min < 0 else wind_range_min
            wind_range_max = (wind_range_max % 360) if wind_range_max > 360 else wind_range_max
            print('HRRR wind angle range: ', (wind_range_min, wind_range_max))
            print('plume wind angle:      ', angle)
    
    else:
        wind_range_min = np.min(wind_array) + degree_incre
        wind_range_max = np.max(wind_array) - degree_incre

        if (0 <= angle <= wind_range_min) or (wind_range_max <= angle <= 360):
            print('Pass wind direction filter.')
            print('HRRR wind angle range: ', (wind_range_min, wind_range_max))
            print('plume wind angle:      ', angle)
            direction = True
        else:
            print('Fail wind direction filter')
            print('HRRR wind angle range: ', (wind_range_min, wind_range_max))
            print('plume wind angle:      ', angle)
        
    return direction


def wind_filter(preserved_file, plume_csv_file, polygon_jsonfile, point_jsonfile, wind_file):
    # Read GeoTIFF image
    preserved_dataset = gdal.Open(preserved_file, gdal.GA_ReadOnly)
    geotransform = preserved_dataset.GetGeoTransform()
    preserved_image = cv2.imread(preserved_file, cv2.IMREAD_UNCHANGED)
    
    # Read csv file
    output_df = pd.read_csv(plume_csv_file)

    # Read GeoJSON data
    with open(point_jsonfile, 'r') as f:
        points_data = json.load(f)
    with open(polygon_jsonfile, 'r') as f:
        polygons_data = json.load(f)
    
    # Extract polygons and points
    polygons = [shape(feature['geometry']) for feature in polygons_data['features']]
    points = [shape(feature['geometry']) for feature in points_data['features']]
    
    # Create output point.json
    filtered_points = []
    wind_direction_list = []

    # Create labels (clumps), label indices should be the same as polygon indices
    preserved_image_binary = ~np.isnan(preserved_image).astype(bool)
    num_labels, labels = cv2.connectedComponents(preserved_image_binary.astype(np.uint8))

    for index in range(1, num_labels):
        polygon = polygons[index - 1]
        print('Prcessing plume index: ', index)

        # Create a sub-image containing a single clump/mask
        _, closest_point = closest_point_to_polygon(polygon, points)
        half_box_size = farthest_dist_point_polygon(polygon, closest_point) + 0.001
        crop_image, start_x, start_y, _, _ = crop_geotiff_image(preserved_image, 
                                                        preserved_file, closest_point.x, closest_point.y, half_box_size)
        crop_labels, _, _, _, _ = crop_geotiff_image(labels,
                                                    preserved_file, closest_point.x, closest_point.y, half_box_size)
        
        xch4a = crop_image.copy()
        xch4a[crop_labels != index] = np.nan
        rows, cols = xch4a.shape

        # Define dat2: 3d array [rows, cols, 3], the 3rd channel: lon, lat, xch4a
        lon_array = np.zeros((rows, cols))
        lat_array = np.zeros((rows, cols))
        for i in range(rows):
            for j in range(cols):
                lon, lat = gdal.ApplyGeoTransform(geotransform, j+start_x, i+start_y)
                lon_array[i, j] = lon
                lat_array[i, j] = lat
        dat2 = np.stack((lon_array, lat_array, xch4a), axis=-1)

        # Intialize iutput raster
        xch4a_sel = np.where(xch4a >= 0, xch4a, np.nan)

        # Intialize output raster
        xch4_final_big = np.full_like(xch4a, np.nan)

        # Create xy_ch4
        xy_xch4 = dat2[:, :, :2]
        xy_xch4 = xy_xch4.reshape(-1, xy_xch4.shape[-1])

        # Create clumps from the input raster
        vv_xch4 = crop_labels.flatten()
        pp_xch4 = xch4a_sel.flatten()

        l_selx = (vv_xch4 == index)
        xch4_final_big = xch4_final_big.flatten()
        xch4_final_big[l_selx] = pp_xch4[l_selx]
        
        # Find center of mass/moment of inertia of XCH4 plume
        cosfac = abs(np.cos(np.mean(xy_xch4[:, 1]) * np.pi / 180))
        X = xy_xch4[:, 0] * cosfac
        Y = xy_xch4[:, 1]
        X0_xch4 = np.mean(np.ma.masked_array(dat2[:, :, 0], mask=np.isnan(xch4a)))
        Y0_xch4 = np.mean(np.ma.masked_array(dat2[:, :, 1], mask=np.isnan(xch4a)))
        xy_centroid = (X0_xch4, Y0_xch4)

        Iyy = np.sum(xch4_final_big[~np.isnan(xch4_final_big)] * (X[~np.isnan(xch4_final_big)] - X0_xch4 * cosfac) ** 2)
        Ixx = np.sum(xch4_final_big[~np.isnan(xch4_final_big)] * (Y[~np.isnan(xch4_final_big)] - Y0_xch4) ** 2)
        Ixy = -np.sum(xch4_final_big[~np.isnan(xch4_final_big)] * (X[~np.isnan(xch4_final_big)] - X0_xch4 * cosfac) *
                        (Y[~np.isnan(xch4_final_big)] - Y0_xch4))

        Mii = np.linalg.eig(np.array([[Ixx, Ixy], [Ixy, Iyy]]))
        i_p = 1 - int(Mii[0][1] > Mii[0][0])
        angle = 180 * np.arctan2(Mii[1][1, i_p], Mii[1][0, i_p]) / np.pi

        Delta_x = np.ptp(X[~np.isnan(xch4_final_big)]) / cosfac
        Delta_y = Delta_x * np.tan(angle * np.pi / 180)

        # Determine exploring distance (unit: meter, 10m mosaic)
        xch4a_count = np.count_nonzero(~np.isnan(xch4a))
        dist = int((np.sqrt(xch4a_count) * 1.5) // 10) * 100

        # Find two ends of the long axis
        dist_ext = dist
        dist_ext_inc = 10
        dat2 = dat2.reshape(-1, dat2.shape[-1])

        pp0 = np.vstack((destPoint(xy_centroid, np.arctan(Delta_x / Delta_y) * 180 / 
                                   np.pi, np.arange(dist_ext, -1, -1 * dist_ext_inc)),
                        destPoint(xy_centroid, np.arctan(Delta_x / Delta_y) * 180 / 
                                  np.pi + 180, np.arange(10, dist_ext + 1, dist_ext_inc))))

        knn = NearestNeighbors(n_neighbors=1).fit(dat2[:, :2])
        ln_nn = knn.kneighbors(pp0, return_distance=False)

        start_ind, end_ind = 0, 0
        for i in range(len(ln_nn)):
            ind = ln_nn[i][0]
            if not np.isnan(dat2[ind, 2]):
                start_ind = ind
                break
        for j in range(len(ln_nn) - 1, -1, -1):
            ind = ln_nn[j][0]
            if not np.isnan(dat2[ind, 2]):
                end_ind = ind
                break

        start_lon, start_lat = xy_xch4[start_ind, 0], xy_xch4[start_ind, 1]
        end_lon, end_lat = xy_xch4[end_ind, 0], xy_xch4[end_ind, 1]

        x_interval_for_plot = [start_lon, end_lon]
        y_interval_for_plot = [start_lat, end_lat]

        # Plot figures
        # Extents for plots
        min_x, max_y = gdal.ApplyGeoTransform(geotransform, start_x, start_y)
        pixel_width = geotransform[1]
        pixel_height = -1 * geotransform[5]
        max_x = min_x + cols * pixel_width
        min_y = max_y - rows * pixel_height
        extent = (min_x, max_x, min_y, max_y)

        plt.figure(figsize=(7, 7))
        plt.imshow(xch4a, extent=extent, cmap='viridis')
        plt.scatter(X0_xch4, Y0_xch4, color='red')
        plt.scatter(start_lon, start_lat, color='blue')
        plt.scatter(end_lon, end_lat, color='black')

        plt.plot(x_interval_for_plot, y_interval_for_plot, color='black')
        plt.title('No. ' + str(index))

        plt.gca().xaxis.set_major_locator(MaxNLocator(4))
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        plt.gca().yaxis.set_major_locator(MaxNLocator(4))
        plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

        
        # Extract HRRR wind data
        wind_df = pd.read_csv(wind_file)
        wind_array = np.array(wind_df.loc[[index - 1]])[0][1:]
        print('HRRR wind array:       ', wind_array)

        # Wind range: +/- 55 degrees
        degree_incre = 55

        # How many points close to this polygon
        dist_thred = 0.001
        points_list = all_points_to_polygon(polygon, points, dist_thred)
        
        # If only one point close to this polygon, this point is plume origin
        if len(points_list) == 1:
            # Add this point to filtered_points
            center_x, center_y = points_list[0][0], points_list[0][1]
            filtered_points.append((center_x, center_y))

            # Among the two ends, find the one that is closer to plume origin
            dist_start2origin = Point(center_x, center_y).distance(Point(start_lon, start_lat))
            dist_end2origin = Point(center_x, center_y).distance(Point(end_lon, end_lat))
            new_xy_plume2 = 'start' if dist_start2origin < dist_end2origin else 'end'

            # Calculate long-axis angle
            if new_xy_plume2 == 'start':
                x_interval = end_lon - start_lon
                y_interval = end_lat - start_lat
            else:
                x_interval = start_lon - end_lon
                y_interval = start_lat - end_lat

            angle = (180 + math.degrees(math.atan2(x_interval, y_interval))) % 360
            print('Origin is start/end:   ', new_xy_plume2)
            
            # Return if wind direction aligns
            direction = wind_angle_align(wind_array, angle, degree_incre)
        
        # If more than one point are close to this polygon
        elif len(points_list) > 1:
            # First, determine the correct long-axis direction
            angle_start = (180 + math.degrees(math.atan2((end_lon - start_lon), (end_lat - start_lat)))) % 360
            angle_end   = (180 + math.degrees(math.atan2((start_lon - end_lon), (start_lat - end_lat)))) % 360
            
            # If one direction is within the wind array range, use this as the long-axis direction
            # Then the point closest to the axis-derived origin is the plume origin
            if (np.min(wind_array) <= angle_start <= np.max(wind_array)) and (angle_end < np.min(wind_array) 
                                                                            or angle_end > np.max(wind_array)):
                new_xy_plume2 = 'start'
                
                _, closest_point = closest_point_to_point(points_list, Point(start_lon, start_lat))
                center_x, center_y = closest_point.x, closest_point.y
                filtered_points.append((center_x, center_y))

                x_interval = end_lon - start_lon
                y_interval = end_lat - start_lat
                angle = (180 + math.degrees(math.atan2(x_interval, y_interval))) % 360
                print('Origin is start/end:   ', new_xy_plume2)

                direction = wind_angle_align(wind_array, angle, degree_incre)

            elif (np.min(wind_array) <= angle_end <= np.max(wind_array)) and (angle_start < np.min(wind_array) 
                                                                            or angle_start > np.max(wind_array)):
                new_xy_plume2 = 'end'
                _, closest_point = closest_point_to_point(points_list, Point(end_lon, end_lat))
                center_x, center_y = closest_point.x, closest_point.y
                filtered_points.append((center_x, center_y))

                x_interval = start_lon - end_lon
                y_interval = start_lat - end_lat
                angle = (180 + math.degrees(math.atan2(x_interval, y_interval))) % 360
                print('Origin is start/end:   ', new_xy_plume2)

                direction = wind_angle_align(wind_array, angle, degree_incre)
            
            # If both directions are out of the wind array range, use the one closer to the range as the long-axis 
            # direction, then the point closest to the axis-derived origin is the plume origin
            elif (angle_start < np.min(wind_array) or angle_start > 
                np.max(wind_array)) and (angle_end < np.min(wind_array) or angle_end > np.max(wind_array)):
                
                angle_start_diff = min(abs(angle_start - np.min(wind_array)), abs(angle_start - np.max(wind_array)))
                angle_end_diff   = min(abs(angle_end - np.min(wind_array)), abs(angle_end - np.max(wind_array)))

                if angle_start_diff < angle_end_diff:
                    new_xy_plume2 = 'start'
                    _, closest_point = closest_point_to_point(points_list, Point(start_lon, start_lat))
                    center_x, center_y = closest_point.x, closest_point.y
                    filtered_points.append((center_x, center_y))
                    x_interval = end_lon - start_lon
                    y_interval = end_lat - start_lat
                else:
                    new_xy_plume2 = 'end'
                    _, closest_point = closest_point_to_point(points_list, Point(end_lon, end_lat))
                    center_x, center_y = closest_point.x, closest_point.y
                    filtered_points.append((center_x, center_y))
                    x_interval = start_lon - end_lon
                    y_interval = start_lat - end_lat
                
                print('Origin is start/end:   ', new_xy_plume2)
                angle = (180 + math.degrees(math.atan2(x_interval, y_interval))) % 360
                direction = wind_angle_align(wind_array, angle, degree_incre)
            
            # If both directions are within the range, use the point closest to one of the axis ends as the plume
            # origin, and the long-axis direction is determined by from that closest end
            else:
                lowest_dist = float('inf')
                for point_x, point_y in points_list:

                    dist_start2origin = Point(point_x, point_y).distance(Point(start_lon, start_lat))
                    dist_end2origin = Point(point_x, point_y).distance(Point(end_lon, end_lat))
                    lower_dist = min(dist_start2origin, dist_end2origin)

                    if lower_dist < lowest_dist:
                        new_xy_plume2 = 'start' if dist_start2origin < dist_end2origin else 'end'
                        closest_point = Point(point_x, point_y)
                        lowest_dist = lower_dist
                
                center_x, center_y = closest_point.x, closest_point.y
                filtered_points.append((center_x, center_y))
                print('Origin is start/end:   ', new_xy_plume2)

                if new_xy_plume2 == 'start':
                    x_interval = end_lon - start_lon
                    y_interval = end_lat - start_lat
                else:
                    x_interval = start_lon - end_lon
                    y_interval = start_lat - end_lat

                angle = (180 + math.degrees(math.atan2(x_interval, y_interval))) % 360
                direction = wind_angle_align(wind_array, angle, degree_incre)
                
        plt.scatter(center_x, center_y, color='red')
        plt.gca().xaxis.set_major_locator(MaxNLocator(4))
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        plt.gca().yaxis.set_major_locator(MaxNLocator(4))
        plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

        wind_direction_list.append(direction)
    
    # Write filtered points to a GeoJSON file
    filtered_features = []
    points_to_keep = {Point(lon, lat): index+1 for index, (lon, lat) in enumerate(filtered_points)}
    for feature in points_data['features']:
        point = Point(feature['geometry']['coordinates'])
        if point in points_to_keep:
            feature['properties']['polygon_index'] = points_to_keep[point]
            filtered_features.append(feature)
    
    filtered_geojson = {'type': 'FeatureCollection', 'features': filtered_features}
    
    filtered_jsonfile = plume_csv_file.split('.csv')[0] + '.geojson'
    with open(filtered_jsonfile, 'w') as f:
        json.dump(filtered_geojson, f)

    output_df['Wind_dir'] = wind_direction_list
    output_df.to_csv(plume_csv_file, index=False)

    print('Plume csv file updated:     ', plume_csv_file)
    print('Point GeoJSON file created: ', filtered_jsonfile)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run wind filter.")
    parser.add_argument("flight", type=str, help="Flight name.")
    
    args = parser.parse_args()
    flight = args.flight

    preserved_file =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '_preserved_v2.tif'
    plume_csv_file =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered_plumes.csv'
    polygon_jsonfile = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered.geojson'
    point_jsonfile =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + 'plumes.geojson'
    wind_file =        '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + '/wind.csv'
    # preserved_file =   '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_preserved_v2_v8.tif'
    # plume_csv_file =   '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' + flight + '_filtered_plumes_test.csv'
    # polygon_jsonfile = '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' + flight + '_filtered_test.geojson'
    # point_jsonfile =   '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' + 'plumes_test.geojson'
    # wind_file =        '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' + '/wind.csv'

    wind_filter(preserved_file, plume_csv_file, polygon_jsonfile, point_jsonfile, wind_file)






# flight = 'RF08_SLC'

# preserved_file =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '_preserved_v2.tif'
# plume_csv_file =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered_plumes.csv'
# polygon_jsonfile = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + flight + '_filtered.geojson'
# point_jsonfile =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + 'plumes.geojson'
# wind_file =        '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + '/wind.csv'

# wind_filter(preserved_file, plume_csv_file, polygon_jsonfile, point_jsonfile, wind_file)

# %%
