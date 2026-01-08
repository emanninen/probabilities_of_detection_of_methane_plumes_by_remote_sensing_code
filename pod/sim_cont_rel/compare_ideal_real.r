
##################
### rough comparison of plumes simulated with idealized vs real wrf les 
##################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(ncdf4)
library(terra)
library(ggplot2)

### read in ideal data
ideal_nc <- nc_open(wrf_les_ideal_file)
real_nc <- nc_open(wrf_les_real_file)

names(ideal_nc$var)
names(real_nc$var)
### coords
long_ideal <- ncvar_get(ideal_nc, "XLONG_V")
long_real <- ncvar_get(real_nc, "XLONG")

lati_ideal <- ncvar_get(ideal_nc, "XLAT_V")
lati_real <- ncvar_get(real_nc, "XLAT")

xch4_ideal <- ncvar_get(ideal_nc, "xch4_1")
xch4_real <- ncvar_get(real_nc, "xch4_1")

### filter missing vals to na
xch4_ideal[xch4_ideal>1e30] <- NA
xch4_real[xch4_real>1e30] <- NA

xch4_mean_ideal <- apply(xch4_ideal, c(1,2),mean, na.rm = T)
# image(xch4_mean_ideal)
xch4_mean_real <- apply(xch4_real, c(1,2),mean, na.rm = T)
# image(xch4_mean_real)

### TODO: recenter plumes
mean_ideal <- rast(xch4_mean_ideal)
### TODO: the ideal wrf les outputs dont have proper lat long data
# ext(mean_ideal) <- c(min(long_ideal), max(long_ideal), min(lati_ideal), max(lati_ideal))
mean_real <- rast(xch4_mean_real)
### recenter plumes to compare
real_center <- xyFromCell( mean_real,which.max(values(mean_real)))
# ext(mean_real)$xmax - real_center[1]
# real_center[2]
mean_real <- crop(mean_real
	, ext(
		real_center[1] - 40.5
		, real_center[1] + 40.5
		, real_center[2] - 40.5
		, real_center[2] + 40.5
))#edn crop
ideal_center <- xyFromCell( mean_ideal,which.max(values(mean_ideal)))
mean_ideal <- crop(mean_ideal
	, ext(
		ideal_center[1] - 40.5
		, ideal_center[1] + 40.5
		, ideal_center[2] - 40.5
		, ideal_center[2] + 40.5
))#edn crop

ext(mean_ideal) <- c(1,dim(mean_ideal)[1],1,dim(mean_ideal)[2])
ext(mean_real) <- c(1,dim(mean_real)[1],1,dim(mean_real)[2])
# plot(mean_real,mean_ideal)

### normalize so can compare sim pluems w diff emiss rates
values(mean_ideal) <- values(mean_ideal)/sum(values(mean_ideal))
values(mean_real) <- values(mean_real)/sum(values(mean_real))

plot_plume(mean_ideal) + labs(title = "Ideal")
plot_plume(mean_real) + labs(title = "Real")
mean_diff <- mean_ideal - mean_real
plot_plume(mean_diff,3)+ labs(title = "Ideal - Real")

compare_df <- data.frame(real = values(mean_real), ideal = values(mean_ideal))
### idk why names wont stick
colnames(compare_df) <- c("real", "ideal")

plot_compare <- ggplot(compare_df) +
	geom_point(aes(x = real, y = ideal), size = 1) +
	geom_abline(slope = 1, intercept = 0) +
	xlim(0,.006) +
	ylim(0,.006) +
	labs(
	     	title = 
		paste(
		      "Real vs Idealized Normalized XCH4, R2:"
		   , round(cor(compare_df$real, compare_df$ideal)^2,3)
		   )
	) +
	ggplot2::theme_classic() +
	ggplot2::theme(text = element_text(size = 30))




# ext(mean_real) <- c(min(long_real), max(long_real), min(lati_real), max(lati_real))

real_center[2]




plot(mean_real,mean_ideal)



plot(as.vector(xch4_mean_real),as.vector(xch4_mean_ideal))

