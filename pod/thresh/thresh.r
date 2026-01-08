
###########################
### for a given MAIR RF
### Sample the variables relevant to PoD (pixel area, gas concentration noise, windspeed)
### So that you can create a visibiility map
### set up to run w/or w/o slurm: sbatch_thresh.sh
###########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(terra)
### params
source(file.path(scripts_dir,"pod/thresh/params.r"))


mos_idx <- 4
mosaic_files <- list.files( list.files(mair_mosaic_dir,full.names = T)[mos_idx] , recursive = T, full.names = T)
mosaic_files <- list.files( mair_mosaic_dir , recursive = T, full.names = T)

mosaic_files <- grep("priority-map", mosaic_files, value = T)
# mosaic_files <- grep("10m", mosaic_files, value = T)
mosaic_files <- grep("nc", mosaic_files, value = T)
### "piceance" has "nc" in it...
mosaic_files <- mosaic_files[!grepl("xml", mosaic_files)]


if(test){mosaic_files <- mosaic_files[1:2]}

mosaic_info <- as.data.frame(do.call(rbind, lapply(1:length(mosaic_files), function(mfi){
	mosaic_file <- mosaic_files[mfi]
	print(mosaic_file)

	info <- scene_info(mosaic_file)
	info$scene_area <- sum(values(cellSize(rast(mosaic_file))))

	vars <- c("area_pixel"
		  , "gas_conc_noise"
		  , "prior_wspd"
		  , "sd_wspd"
		  , "scene_file"
		  , "time_start"
		  , "time_end"
		  , "utc"
		  , "scene_area"
		  , "p_surf"
	)#ednd c()
	return(as.data.frame(info[vars]))
})))#end .alpply rbind

write.csv(
	mosaic_info
	, ifelse(
		 test
		 , file.path(thresh_dir, "thresh_test.csv")
		 , file.path(thresh_dir, "thresh.csv")
	)#edn ifelse
	, row.names = F 
)#end wrt
