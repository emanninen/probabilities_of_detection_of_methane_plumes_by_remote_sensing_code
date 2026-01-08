
###########################
### get characterist gas conc noise & pixel area MSAT
###########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(terra)
### params
source(file.path(scripts_dir,"pod/thresh/params.r"))

### msat files
msat_files <- list.files(msat_ss_dir , recursive = T, full.names = T)

msat <- rast(msat_files[1])

msat_scene_info <- scene_info(msat_files[1])

msat_scene_info$time_start
msat_scene_info$area_pixel
msat_scene_info$gas_conc_noise
msat_scene_info$prior_wspd

