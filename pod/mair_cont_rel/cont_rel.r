
#########################
### take MAIR controlled release scenes and extract SNR vars + format them to be beta-binomiales
### input
### ### list of filenames of controlled release scenes
### output
### ### saved rds files containing all the info needed to run beta-binomial on the SNR
### depends
### ### trial_scene()
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(terra)
library(ncdf4)
library(lubridate)
### params
source(file.path(scripts_dir,"pod/mair_cont_rel/params.r"))

##########################
### natty convert wrf les plumes
##########################

### clear trial dir
file.remove(list.files(mair_cont_rel_trial_natty_dir, full.names = T))

for(rf_idx in 1:length(rfs)){
	rf <- rfs[rf_idx]
	list.files(mair_scenes_dir)
	cont_rel_files <- list.files(mair_scenes_dir, full.names = T, pattern = rf)
	###  pic release location
	cont_rel_location <- cont_rel_locations[[rf]]

	for(seg_idx in 1:length(cont_rel_files)){

		mair_file <- cont_rel_files[seg_idx]

		trial <- trial_scene(
			scene_file = mair_file
			, cont_rel_location = cont_rel_location
			, meter_file = file.path(mair_cont_rel_meter_emiss_rate_dir, "ju_meter.csv")
			, rf = rf
			, trial_dir = mair_cont_rel_trial_natty_dir
			, save = T
		)#end mair_cont_rel_segment
	}#end if F
}#end for

### dont enhance natty non detects
### run di detection on the cont rel

### PHEW bullshit of moving the cont rel files
trial_files <- list.files(mair_cont_rel_trial_natty_dir, full.names = T)
# readRDS(trial_files[1])$scene_file
# scene_files <- sapply(trial_files,function(ff){ readRDS(ff)$scene_file })
# new_scene_files <- file.path(mair_scenes_dir,paste0(sapply(trial_files,function(ff){ round(readRDS(ff)$area_pixel)}),str_split_i(scene_files,"/",14)))
# file.remove(list.files(mair_scenes_dir,full.names = T))
# file.copy(scene_files,new_scene_files, overwrite = T)

cont_rel <- di_detection(
	trial_files
	, detection_dir = mair_detection_dir
	, detection_file = "natty_di_detection.csv"
	, save_cont_rel = T
)#end di detect

cont_rel$DD 
natty_no_detect_files <- cont_rel[!cont_rel$DD,"scene_file"]

##########################
### enhanced convert wrf les plumes
##########################

### clear trial dir
file.remove(list.files(mair_cont_rel_trial_enhanced_dir, full.names = T))

for(rf_idx in 1:length(rf_dirs)){

	rf <- rfs[rf_idx]
	print(paste("rf:", rf))
	list.files(mair_scenes_dir)
	cont_rel_files <- list.files(mair_scenes_dir, full.names = T, pattern = rf)

	###  pic release location
	cont_rel_location <- cont_rel_locations[[rf]]

	for(seg_idx in 1:length(cont_rel_files)){

		mair_file <- cont_rel_files[seg_idx]
		### dont try to enhance when dont detect natty
		if(mair_file %in% natty_no_detect_files){next}

		print(mair_file)
		trial <- trial_scene(
			scene_file = mair_file
			, cont_rel_location = cont_rel_location
			, meter_file = file.path(mair_cont_rel_meter_emiss_rate_dir, "ju_meter.csv")
			, rf = rf
			, enhance = T

			, target_emiss = qq_samples
			, target_noise = nn_samples
			, target_area_pixel_fact = apf_samples

			, trial_dir = 
				ifelse(
					test
					, mair_cont_rel_trial_enhanced_test_dir
					, mair_cont_rel_trial_enhanced_dir
				)#end if else
			, save = T
		)#end mair_cont_rel_segment
	}#end if F
}#end for
