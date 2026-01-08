
################################
### examine jus wrf les files to see
### how to make a test images for training beta binom
################################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(ncdf4)

wrf_les_ideal_outputs_dir <- "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/WRF/WRF_outputs/LES_outputs/LES_2024_v2/"
scene_files <- list.files(wrf_les_ideal_outputs_dir, full.names = T)

wrf_file <- scene_files[1]
# wrf_file <- list.files(wrf_les_dir, full.names = T)[2]

wrf_nc <- nc_open(wrf_file)
names(wrf_nc$var)
### looking for pressure var
pb <- ncvar_get(wrf_nc, "PB")
ph <- ncvar_get(wrf_nc, "PH")
phb <- ncvar_get(wrf_nc, "PHB")
psfc <- ncvar_get(wrf_nc, "PSFC")

### lat lon
### there is an entry for each time hopefully lat lon time invariant
lon <- ncvar_get(wrf_nc, "XLONG")[,,1]
lat <- ncvar_get(wrf_nc, "XLAT")[,,1]

### trace is concentration in each grid cell 
# tracer_1 <- ncvar_get(wrf_nc, "tracer_1")
### for now only take xch41
### xch4 2-5 i dont understand
xch4_1 <- ncvar_get(wrf_nc, "xch4_1")

### wind
uu <- ncvar_get(wrf_nc, "U")
vv <- ncvar_get(wrf_nc, "V")

tt <- 50
### wind level TODO: ask what I should take
wind_level <- 2
### I think the plume is centered...
scene_plume_loc <- c(60,60)
### TODO: look at the plume location, 10m winds
wspd <- sqrt(
	uu[scene_plume_loc[1],scene_plume_loc[2],wind_level,tt]^2
	+ vv[scene_plume_loc[1],scene_plume_loc[2],wind_level,tt]^2
)#end squrt

hist(wspd)

input_emiss_rate <- 100
target_emiss_rate <- 500  

emiss_rate_scaling <- target_emiss_rate/input_emiss_rate

xch4 <- xch4_1[,,tt]
image(xch4)
### scale xch4 bc it should be linear w emission rate
xch4_scaled <- xch4*emiss_rate_scaling
image(xch4_scaled )


hist(xch4_scaled)
gas_conc_noise <- 100e-3
### add noise
xch4_noise <- xch4_scaled+rnorm(dim(xch4_scaled)[1]*dim(xch4_scaled)[2], 0, gas_conc_noise)
image(xch4_noise)

cont_rel_out_dir <- "/n/wofsy_lab2/Users/emanninen/pod/sim_cont_rel"

lon <- ncdim_def("lon",units = "deg", lon)
lat <- ncdim_def("lat",units = "deg", lat)
emiss_rate <- ncvar_def("emiss_rate", units = "kghr-1", dim = 1)

trial <- nc_create(
	file.path(cont_rel_out_dir, paste(tt,target_emiss_rate,round(wspd,2), round(gas_conc_noise,2), ".nc",sep = "_") )
	, vars = 
)#end nc_create


