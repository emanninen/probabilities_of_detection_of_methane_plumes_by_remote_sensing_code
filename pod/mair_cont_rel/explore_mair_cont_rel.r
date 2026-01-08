
### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")

### dpeends
library(terra)
library(ncdf4)
library(lubridate)

### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)

rf_dirs <- list(mair_cont_rel_rf04, mair_cont_rel_rf05)
cont_rel_dirs <- list.files(mair_cont_rel_rf04, full.names = T)[3:13]
cont_rel_files <- list.files(cont_rel_dirs, full.names = T, recursive = T)
cont_rel_files <- grep(".nc", cont_rel_files, value = T)

ten_files <- grep("10m", cont_rel_files, value = T)
thirty_files <- grep("30m", cont_rel_files, value = T)

cont_rel_locations <- list( 
	rf04 <- c(-102.29, 32.056)
	, rf05 <- c(-102.29, 32.056)
)#end list

### clear trial dir
file.remove(list.files(mair_cont_rel_trial_dir, full.names = T))

rf_idx <- 1
seg_idx <- 1
for(rf_idx in 1:2){
	print(rf_idx)
	rf_dir <- rf_dirs[[rf_idx]]
	cont_rel_files <- list.files(rf_dir, full.names = T)
	# ten_files <- grep("10m", cont_rel_files, value = T)
	# thirty_files <- grep("30m", cont_rel_files, value = T)
	mair_files <- c(ten_files, thirty_files) 
	for(seg_idx in 1:length(mair_files)){
		mair_file <- mair_files[seg_idx]
		print(mair_file)
	}#end forj
}#end for

