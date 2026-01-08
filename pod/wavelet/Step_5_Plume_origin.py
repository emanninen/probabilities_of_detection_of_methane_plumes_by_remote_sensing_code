import os
import argparse
import numpy as np
from netCDF4 import Dataset
from math import atan2
from pyproj import Geod
from skimage import measure
from skimage.measure import inertia_tensor, regionprops
from skimage.measure._regionprops import RegionProperties
from msat_level4_divergence_integral.io.hrrr_input import extract_met_data, extract_met_data_from_nc
from msat_level4_divergence_integral.io.layer_averaging import select_layer
from msat_level4_divergence_integral.divergence_integral.integration import prep_met, _prep_ch4
from msat_level4_divergence_integral.io.l3_input import format_level3
from msat_level4_divergence_integral.plume_finder.plume_characterization import *
from msat_level4_divergence_integral.plume_finder.types import PlumeRotationSource, PlumeFinderOptions
from msat_level4_divergence_integral.divergence_integral.types import DivergenceIntegralOptions
from msat_level4_divergence_integral.plume_pruning.plume_pruning import prune_plumes
from msat_level4_divergence_integral.io.geojson import write_plumes
geod = Geod(ellps="WGS84")


def _weighted_orientation(clump: RegionProperties) -> float:
    mu = clump.moments_weighted_central
    tensor = inertia_tensor(np.where(clump.image > 0, clump.image_intensity, 0), mu=mu)
    
    # Formula A
    # Copied from `skimage/measure/_regionprops.py`
    a, b, b, c = tensor.flat

    if a - c == 0:
        if b < 0:
            return np.pi / 4.0
        else:
            return -np.pi / 4.0
    else:
        return 0.5 * atan2(-2 * b, c - a)
    
def _clump_derived_wind_angle(clump: RegionProperties, met_wind_angle: float) -> float:
    # This will take on a value between 0˚ - 180˚
    derived_angle = np.rad2deg(np.mod(np.pi / 2 - _weighted_orientation(clump), 2 * np.pi))

    # Find the smallest angle between the plume-derived angle and the wind angle
    # given by the met model
    derived_vs_met = min(
        np.abs(met_wind_angle - derived_angle), 360 - np.abs(met_wind_angle - derived_angle)
    )

    # Rotate the plume-derived angle by 180˚ if the two angles differ by >90˚
    if derived_vs_met > 90:
        derived_angle = np.mod(derived_angle + 180, 360)

    return derived_angle

def _wind_angle(coord: np.ndarray, wind_interpolator: LinearNDInterpolator) -> float:
    # east wind, north wind
    u, v = np.split(wind_interpolator(coord[np.newaxis, :]), 2, axis=-1)

    # arctan2 signature is y,x -> theta
    # rotate by pi (transition from "blowing towards" to "blowing from"), then convert to degrees
    thetas_rad = np.mod(-1 * np.arctan2(v, u) + 3 * np.pi / 2, 2 * np.pi)
    thetas_deg = np.rad2deg(thetas_rad)

    return thetas_deg[0]

def index_2_lonlat(ch4_vertices, x, y):
     return ch4_vertices[x, y, :][0].item(), ch4_vertices[x, y, :][1].item()

def determine_plume_origin(clump: RegionProperties, wind_interpolator: LinearNDInterpolator, ch4_vertices):
    # Extract clump eccentricity
    eccentricity = clump.eccentricity
    
    # Extract clump centroid coordinates
    centroid_x, centroid_y = clump.centroid
    lon, lat = index_2_lonlat(ch4_vertices, round(centroid_x), round(centroid_y))
    
    # Extract clump pixel coordinates
    coords = clump.coords
    clump_coords = np.array([index_2_lonlat(ch4_vertices, round(x), round(y)) for x, y in coords])
    
    # Extract HRRR wind angle based on centroid coordinates
    met_wind_angle = _wind_angle(np.array((lon, lat)), wind_interpolator)[0]

    # Create clump-derived wind angle for the clumps with high eccentricity
    if eccentricity > PlumeFinderOptions.xch4_wind_clump_eccentricity_threshold:
        plume_wind_angle = _clump_derived_wind_angle(clump, met_wind_angle)
    else:
        plume_wind_angle = None
    
    # Determine wind angle
    wind_angle = met_wind_angle if plume_wind_angle is None else plume_wind_angle
    
    # Determine plume origin and end
    plume_origin = furthest_upwind_location(clump_coords, wind_angle)
    plume_end = furthest_upwind_location(clump_coords, (wind_angle + 180) % 360)
    
    if plume_wind_angle is not None:
        rotation_angle_source = PlumeRotationSource.plume
    else:
        rotation_angle_source = PlumeRotationSource.weather_model
    
    return Plume(
        x=plume_origin[0],
        y=plume_origin[1],
        elev=None,  # Not using right now
        wind_speed=None,  # Not using right now
        wind_angle=met_wind_angle,
        plume_eccentricity=eccentricity,
        plume_length=geod.line_length(
            lons=[plume_origin[0], plume_end[0]], lats=[plume_origin[1], plume_end[1]]
        ),
        plume_wind_angle=plume_wind_angle,
        high_xch4=len(clump_coords) > 0,
        rotation_angle_source=rotation_angle_source,
    )

def plume_origin(directory_path, output_path, l3_30m_file, hrrr_file):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    # Prepare wind data
    hrrr_path = os.path.join(directory_path, hrrr_file)
    met_data = extract_met_data(hrrr_path)
    # met_data = extract_met_data_from_nc(hrrr_path)
    layer_met_data = select_layer(met_data)
    prepared_met_data = prep_met(layer_met_data)

    wind_interpolator = prepared_met_data.wind_interpolator

    # Prepare xch4 data
    l3_30m_path = os.path.join(directory_path, l3_30m_file)
    l3_30m_dataset = Dataset(l3_30m_path, "r")
    l3_xch4_30m = format_level3(l3_30m_dataset, num_samples_threshold=
                                    DivergenceIntegralOptions.num_samples_threshold)
    ch4_vertices, ch4 = _prep_ch4(l3_xch4_30m)
    
    # Create clumps based on connected component algorithm
    ch4_image = np.squeeze(ch4, axis=-1)
    binary_image = ~np.isnan(ch4_image)
    labeled_image = measure.label(binary_image, connectivity=2)
    properties = measure.regionprops(labeled_image, intensity_image = l3_xch4_30m[:,:,2])

    # Determine plume origins for all the clumps
    plumes = [determine_plume_origin(clump, wind_interpolator, ch4_vertices) for clump in properties]

    if len(plumes) > 0:
        print(f"Plume finder complete, {len(plumes)} plumes found for ", flight)
    else:
        print("Plume finder complete, No plumes were found for ", flight)

    # # Prune plumes
    # in_pruned_set_idxs = prune_plumes(plumes, properties, l3_xch4_30m)
    # print("Pruning complete")
    in_pruned_set_idxs = list(range(len(plumes)))

    # Write plumes and clumps to output directory
    write_plumes(plumes, output_path + "/plumes.geojson", in_pruned_set_idxs)






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine plume origin.")
    parser.add_argument("flight", type=str, help="Flight name.")
    parser.add_argument("hour", type=str, help="Hour of HRRR data.")
    parser.add_argument("file_name", type=str, help="Masked 30m NetCDF file name.")
    
    args = parser.parse_args()
    flight = args.flight
    hour = args.hour
    l3_30m_file = args.file_name
    
    hrrr_file =        'hrrr.t' + hour + 'z.wrfnatf.grib2'
    directory_path =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/'
    output_path =      '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/'
    # hrrr_file =        'hrrr.t' + hour + 'z.wrfnatf.grib2'
    # directory_path =   '/n/home12/zhanzhang/Test_DI/' + flight + '/'
    # output_path =      '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/'

    plume_origin(directory_path, output_path, l3_30m_file, hrrr_file)



# flight =         'MSAT_240813'
# hrrr_file =      '20240522_18-23_hrrr_21.nc'
# l3_30m_file =    'MethaneSAT_L3_45m_20240522T211158_20240522T211230_chairfield_upd_v7.nc'
# directory_path = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/'
# output_path =    '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_v2/'

# plume_origin(directory_path, output_path, l3_30m_file, hrrr_file)