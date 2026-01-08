#%%
import os, sys
import shutil
from netCDF4 import Dataset
from Image_proc_functions import *
import argparse
import matplotlib.pyplot as plt



# Create a new NetCDF 30m.nc file based on the masked image
# Input nc_file: 10m.nc
def update_nc_input(masked_file, nc_file):
    # Read masked GeoTIFF image
    masked_dataset = gdal.Open(masked_file, gdal.GA_ReadOnly)
    masked_image = masked_dataset.GetRasterBand(1).ReadAsArray()

    # Shrink 10m GeoTIFF image to 30m
    masked_image_30m = shrink_image(masked_image)
    # masked_image_30m = masked_image
    masked_image_30m = np.flip(masked_image_30m, axis=0)
    
    # Copy the 30m.nc to be the updated output
    nc_file_30m = nc_file.split('_10m.nc')[0] + '_30m.nc'
    # nc_file_30m = nc_file
    upd_nc_file_30m = nc_file_30m.split('.nc')[0] + '_upd.nc'
    shutil.copyfile(nc_file_30m, upd_nc_file_30m)

    upd_nc_30m = Dataset(upd_nc_file_30m, 'r+')
    upd_xch4_30m = upd_nc_30m.variables['xch4'][:]
    
    # Set pixels masked out in the masked image as NAN in the updated output
    binary_mask = np.isnan(masked_image_30m).astype(bool)
    upd_xch4_30m_masked = np.ma.MaskedArray(upd_xch4_30m, mask=binary_mask)

    upd_nc_30m.variables['xch4'][:] = upd_xch4_30m_masked
    upd_nc_30m.close()
    
    print("NetCDF 30m file updated:", upd_nc_file_30m)






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update NetCDF Inputs (30m.nc).")
    parser.add_argument("flight", type=str, help="Flight name.")
    parser.add_argument("file_name", type=str, help="10m NetCDF file name.")
    
    # Run the function with the provided arguments
    args = parser.parse_args()
    flight = args.flight
    file_name = args.file_name

    masked_file =    '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/Output_wavelet/' + flight + '_preserved.tif'
    directory_path = '/n/holylfs04/LABS/wofsy_lab/Lab/zzhang/' + flight + '/'
    # masked_file =    '/n/home12/zhanzhang/Test_DI/' + flight + '/Output_wavelet/' + flight + '_preserved_v2_v7.tif'
    # directory_path = '/n/home12/zhanzhang/Test_DI/' + flight + '/'

    nc_file = os.path.join(directory_path, file_name)

    update_nc_input(masked_file, nc_file)