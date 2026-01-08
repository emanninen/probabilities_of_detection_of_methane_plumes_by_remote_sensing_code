# Plume masking using wavelet denoising

## Overview
This repository contains the code for point source plume masking based on the denoised plume concentration image using wavelet transform techniques. Given a plume concentration image, the code extracts the high-frequency pixels by conducting wavelet transform and inverse wavelet transform. Then it subtracts the high-frequencies from the original image, assuming that the high-frequencies are mainly background noise. The resulting masks then go through three main filters: "hotspot" filter, wind direction filter, and shape filter. The masks passing all the filters are preserved in the end. The code also determines plume origin locations and quantifies emission rates using the plume finding algorithm and growing box algorithm from the Divergence Integral method.

## Contents

1. **Step 1: Run wavelet**

This step denoises a plume concentration image by extracting and subtracting the high-frequency pixels using wavelet transform techniques. First, if the original image has NAN-value holes, the code fills the holes with random Gaussian noise to have a complete image. Next, the code sets strong signals (mu + scale * sigma) to consistent values (max pixel value), so that these strong signals are less likely to be discarded as high-frequencies. Then, the high-frequency pixels are extracted and subtracted by running wavelet transform and inverse wavelet transform. Further wavelet denoising based on soft thresholding is also conducted. The resulting image is the denoised version of the original image. Finally, based on the pixel value threshold (mu + sigma_factor * sigma) and the connected component algorithm, plume masks are generated from the denoised image (added with Gaussian filtering).

- Code files: Step_1_Run_wavelet.py, Wavelet_functions.py
- Input: flight_name, NetCDF plume concentration image
- Output: 
  - Original GeoTIFF image: flight_name.tif
  - Input GeoTIFF image (filled NAN holes): flight_name_filled.tif
  - Reconstruncted image (high frequencies)
  - Subtract image (input - reconstructed)
  - Denoised GeoTIFF image (further denoised): flight_name_denoised.tif
  - Masked GeoTIFF image: flight_name_masked.tif
- Parameters:
  - scale: factor of the consistent signal threshold
  - sigma_factor: factor of the plume masking threshold

2. **Step 2: Mask processing (regional wavelet)**

This step aims to do regional wavelet processing to split the large masks created in the previous step into smaller ones. Among all the clumps masked by the previous step, the ones with size greater than the large_label_mask_threshold are considered as input, and go through another round of wavelet processing. Note that a resulted binary mask may have False values within it, which leads to NAN values within a potential plume clump. To address this, this step also fills the NAN-value holes within the clumps with values from the original image.

- Code files: Step_2_Mask_processing.py
- Input: flight_name, Original GeoTIFF image: flight_name.tif, Masked GeoTIFF image: flight_name_masked.tif
- Output: 
  - A new version of masked GeoTIFF image: flight_name_masked_v2.tif
- Parameters:
  - scale: factor of the consistent signal threshold
  - sigma_factor: factor of the plume masking threshold
  - large_label_mask_threshold: the maximum mask size threshold

3. **Step 3: Hotspot filtering**

This step runs "Hotspot" filter to the masks created in the previous step to filter out those with low concentrations. "Hotspot" is defined as a smaller clump (size greater than hotspot_label_mask_threshold) with strong signals (value higher than hotspot_threshold) within the mask. Only the masks with valid hotspots within them are preserved after this step. The value of hotspot_threshold is clump specific, which is calculated based on the standard deviation of the background surrounding the clump. First, from a given percentile threshold (base_hotspot_percentile, default 99.97%) and the pixel values of the whole flight scene, we calculate the base hotspot_threshold (base_hotspot_threshold) and the background standard deviation (base_background_std). Next, we define a sub-region surrounding the clump and calculate its background standard deviation (sub_background_std). Then, if (sub_background_std - base_background_std) > 2, the sub-region is considered as noisy, so we want to raise hotspot_threshold to avoid picking noise signal as plume signal. The hotspot_threshold for this clump thus equals to base_hotspot_threshold + 3.4 * (sub_background_std - base_background_std), as (mu + 3.4 * sigma) approximately equals to 99.97% in the Gaussian distribution. In the end, this code also filters out small masks (size smaller than label_mask_threshold).

- Code files: Step_3_Hotspot.py
- Input: flight_name, Original GeoTIFF image: flight_name.tif, Masked GeoTIFF image: flight_name_masked_v2.tif
- Output: 
  - A new version of masked GeoTIFF image filtered by base_hotspot_threshold: flight_name_masked_v3.tif
  - A new version of masked GeoTIFF image filtered by hotspot_threshold (updated based on sub-region std): flight_name_preserved.tif
- Parameters:
  - label_mask_threshold: the minimum mask size threshold
  - hotspot_label_mask_threshold: the minimum hotspot size threshold
  - base_hotspot_percentile: the percentile threshold of the whole flight scene

4. **Step 4: Update input**

This step converts the GeoTIFF image into the NetCDF format, in order to be used in the next step.

- Code files: Step_4_UpdateNCInput.py
- Input: flight_name, NetCDF file: xxx_10m.nc, Masked GeoTIFF image: flight_name_preserved.tif
- Output: NetCDF file of the masked flight scene: xxx_30m_upd.nc

5. **Step 5: Determine plume origin**

This step runs plume finder code from L4-DI method to determine plume origins.

- Code files: Step_5_Plume_origin.py
- Input: flight_name, NetCDF file of the masked flight scene: xxx_30m_upd.nc, HRRR wind file: xxx.grib2
- Output: GeoJSON file of plume origin points: plumes.geojson

6. **Step 6: Create initial mask results**

This step preserves the masks with plume origins nearby into a GeoJSON file. Since the plume origins created in the previous step are usually either at the mask boundary or very close to the masks, usually all the masks pass this step. The GeoJSON file also contains attribute tables with plume index, lat, lon, hotspot_threshold, sub_background_std, mask size, hotspot size, hotspot size ratio (hotspot size / mask size), fibre length ratio, and fibre width. The latter two parameters are the main parameters of the shape filter. They aim to discard the masks with the "spider" shape, i.e., having multiple long thin branches spreading out. To do so, we calcualte both the long-axis length and fibre length of each mask. Then, we define the fibre length ratio as fibre length / long-axis length. The higher the fibre length ratio is, the more likely that the mask has a "spider" shape. In the meantime, we define the fibre width as area / fibre length. The lower the fibre width is, the more likely that the mask has a "spider shape". 

- Code files: Step_6_Post_processing.py
- Input: flight_name, Original GeoTIFF image: flight_name.tif, Masked GeoTIFF image: flight_name_preserved.tif, GeoJSON file of plume origins: plumes.geojson
- Output: 
  - GeoJSon file of plume polygons: flight_name_filtered.geojson
  - CSV file of plumes: flight_name_filtered_plumes.csv
- Parameters:
  - base_hotspot_percentile: the percentile threshold of the whole flight scene

7. **Step 7: Extract HRRR wind data**

This step extracts HRRR wind angles for all the masks preserved in the previous step, in order to be used in the wind filter step. The functions are modified from the L4-DI code.

- Code files: Step_7_Wind.py
- Input: flight_name, GeoJSON file of plume polygons: flight_name_filtered.geojson
- Output: CSV file of HRRR wind: wind.csv

8. **Step 8: Quantify emission rates**

This step runs growing box code from L4-DI method to quantify the emission rates for all the plume origins.

- Code files: Step_8_Quantify.py
- Input: flight_name, GeoJSON file of plume origins: plumes.geojson
- Output: GeoJSON file of emission rates: growing_boxes.geojson

9. **Step 9: Wind filter**

This step runs wind direction filter to the masks to filter out those with main-axis direction contradicting to the HRRR wind direction. The code extracts HRRR wind angles for each mask from Step 7, and then form an angle range by adding angle buffers (+/-55). Then, only the masks with main-axis direction within the HRRR wind range are preserved. The code also chooses the plume origin if multiple origins are determined for the same mask in Step 5. In this case, the point closest to the upwind end of the plume is chosen as the plume origin. The code adds a binary column of "wind_filter" to the plume csv file, and creates a new version of plume origin GeoJSON file.

- Code files: Step_9_Wind_filter.py
- Input: 
- Output: GeoJSon file of filtered plume origins: flight_name_filtered_plumes.geojson

10. **Step 10: Final filter**

This step does final mask filtering based on the paramter values of "hotspot" filter and shape filter. It also filters out the maks with too many NAN pixels along the boundary. This boundary_nan_filtering method aims to discard those masks created along the boundary of NAN-value holes in the original image. These NAN-value holes can be, for example, cloud and cloud shadows.

- Code files: Step_10_Final_filter.py
- Input: 
- Output: 
  - A new version of plume polygn GeoJSon file: flight_name_filtered_v2.geojson
  - A new version of plume origin GeoJSon file: flight_name_filtered_plumes_v2.geojson
  - A new version of plume CSV file: flight_name_filtered_plumes_v2.csv
