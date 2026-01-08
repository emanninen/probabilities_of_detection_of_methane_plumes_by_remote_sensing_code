import os
import argparse
from netCDF4 import Dataset
from msat_level4_divergence_integral.io.geojson import (
    read_plumes,
    write_growing_boxes,
    write_growing_box_details,
)
from msat_level4_divergence_integral.divergence_integral.integration import (
    prep_met,
    compute_divergence_integral_growing_box,
)
from msat_level4_divergence_integral.io.l3_input import format_level3
from msat_level4_divergence_integral.io.hrrr_input import extract_met_data, extract_met_data_from_nc
from msat_level4_divergence_integral.io.layer_averaging import select_layer
from msat_level4_divergence_integral.growing_box.box_building import gen_growing_box_sets
from msat_level4_divergence_integral.divergence_integral.types import DivergenceIntegralOptions
from msat_level4_divergence_integral.growing_box.types import GrowingBoxOptions


def quantify_plumes(directory_path, output_path, l3_10m_file, hrrr_file, point_jsonfile):
    # Prepare wind data
    hrrr_path = os.path.join(directory_path, hrrr_file)
    met_data = extract_met_data(hrrr_path)
    # met_data = extract_met_data_from_nc(hrrr_path)
    layer_met_data = select_layer(met_data)
    prepared_met_data = prep_met(layer_met_data)

    # Prepare xch4 data
    l3_10m_path = os.path.join(directory_path, l3_10m_file)
    l3_10m_dataset = Dataset(l3_10m_path, "r")
    l3_xch4_10m = format_level3(l3_10m_dataset, num_samples_threshold=
                                    DivergenceIntegralOptions.num_samples_threshold)

    # Read point geojson file
    sorted_plumes, in_pruned_set_plume_ids = read_plumes(point_jsonfile)

    # Calculate emission rates
    gb_options = GrowingBoxOptions()
    growing_box_sets = gen_growing_box_sets(sorted_plumes, gb_options)

    growing_box_results = compute_divergence_integral_growing_box(
            growing_box_sets=growing_box_sets,
            ch4_data=l3_xch4_10m,
            prepared_met_data=prepared_met_data,
            options=gb_options,
        )

    # Write growing box output to output path
    write_growing_boxes(
        growing_box_results = growing_box_results,
        growing_box_sets = growing_box_sets,
        local_path = output_path + "/growing_boxes.geojson",
        in_pruned_set_idxs = in_pruned_set_plume_ids,
    )
    # write_growing_box_details(
    #     growing_box_results = growing_box_results,
    #     growing_box_sets = growing_box_sets,
    #     local_path = output_path + "/growing_box_details.geojson",
    # )
    
    print('Quantification done. Output file: ', output_path + "/growing_boxes.geojson")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quantify emission rates.")
    parser.add_argument("flight", type=str, help="Flight name.")
    parser.add_argument("hour", type=str, help="Hour of HRRR data.")
    parser.add_argument("file_name", type=str, help="10m NetCDF file name.")
    
    args = parser.parse_args()
    flight = args.flight
    hour = args.hour
    l3_10m_file = args.file_name
    
    hrrr_file =        'hrrr.t' + hour + 'z.wrfnatf.grib2'
    directory_path =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/'
    point_jsonfile =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/' + 'plumes.geojson'
    output_path =      '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output/'
    # hrrr_file =        'hrrr.t' + hour + 'z.wrfnatf.grib2'
    # directory_path =   '/n/home12/zhanzhang/Test_DI/' + flight + '/'
    # point_jsonfile =   '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/' + 'plumes_test.geojson'
    # output_path =      '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_v2/'

    quantify_plumes(directory_path, output_path, l3_10m_file, hrrr_file, point_jsonfile)




# flight =           'MSAT_240813'
# l3_10m_file =      'MethaneSAT_L3_45m_20240522T211158_20240522T211230_chairfield.nc'

# hrrr_file =        '20240522_18-23_hrrr_21.nc'
# directory_path =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/'
# point_jsonfile =   '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_v2/' + 'plumes_test.geojson'
# output_path =      '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_v2/'

# quantify_plumes(directory_path, output_path, l3_10m_file, hrrr_file, point_jsonfile)