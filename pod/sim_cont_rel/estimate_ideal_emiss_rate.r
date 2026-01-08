
##################################
### try to estimate the emissions rate of idealized wrf
##################################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)


scene_files <- list.files(wrf_les_ideal_outputs_dir, full.names = T)
scene_file <- scene_files[8]
print(scene_file)

### load in data nc
wrf_nc <- ncdf4::nc_open(scene_file)
### load in xch4 data full time series 
names(wrf_nc$var)
xch4_varname <- "tr17_1" 
xch4_varname <- "xch4_1" 
xch4_full <- ncdf4::ncvar_get(wrf_nc, xch4_varname)

xch4_T <- ncdf4::ncvar_get(wrf_nc, "T")
xch4_time <- ncdf4::ncvar_get(wrf_nc, "Times")
### Ju didn't sum the vertical dim for this batch
# xch4_full <- apply(xch4_full, c(1,2,4), sum) 
### look at first couple time steps
xch4_0 <- xch4_full[,,1]
image(xch4_0)
sum(xch4_0)
xch4_1 <- xch4_full[,,2]
sum(xch4_1)
image(xch4_1)
xch4_2 <- xch4_full[,,3]
image(xch4_2)
### looks like zero ch4 in the air at first time step, then theres some..
### hopefully we can determine how much got added in each time step

total_xch4 <- vector()
for(ti in  1:length(xch4_time)){
	xch4 <- xch4_full[,,ti]
	total_xch4 <- c(total_xch4, sum(xch4))
}#end for 1:length xc4 time

delta_xch4 <- diff(total_xch4)

hist(delta_xch4)
delta_xch4[1]
### only look at the start so not affected by outflow
plot(2:length(xch4_time),delta_xch4)

plot(2:11,delta_xch4[1:10])

### convert inflow in ppb (?)/scene/min to kg/scene/hr 
inflow_xch4 <- mean(delta_xch4[1:5])

m2_pixel = 33.3^2 # m2

m2_scene <- #m2
	m2_pixel *
	dim(xch4_full)[1] *
	dim(xch4_full)[2] 

min_hr <- 60 # min/hr
kgCH4_molsCH4 <- 16.042e-3 #kgch4 /mol ch4
p_surf <- mean(ncdf4::ncvar_get(wrf_nc, "PSFC")) # Pa
kgdryair_m2 <- kgdryair(p_surf) # 
m2_kgdryair <- 1/kgdryair(p_surf)
p_ppm <- 1e-6
p_ppb <- 1e-9
kgdryair_moldryair <- 28.96e-3
moldryair_kgdryair <- 1/28.96e-3
# 39km

### mass correction: wrf les goes up to 79e3 Pa
### THIS IS THE DIFFERENCE 
column_mass_correction <- (kgdryair(100e3)-kgdryair(79e3))/kgdryair(100e3)

### WINNER !!! 
###
734.93
###
###


# molair_scene <- 
# 	m2_scene *
molair_pixel <- 
	m2_pixel *
	kgdryair_m2 *#kgdryair_scene
	moldryair_kgdryair * #moldryair_scene
	column_mass_correction
	
kgch4_pch4 <- 
	molair_pixel * 
	# molair_scene * 
	kgCH4_molsCH4

ppm_minscene <- inflow_xch4

###
kgch4_hrscene <- ppm_minscene *
	min_hr *# ppm_hrscene
	p_ppm *# p ch4_hrscene
	kgch4_pch4 # kgch4_hrscene 

print(paste("emiss rate kg/hr:",kgch4_hrscene))

