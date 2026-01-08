
#########################
### run wave detection
### take controlled release scenes (mair + wrfles), run the wave detection method 
### and get controlled release data with detect/no detect
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)

### params
source(file.path(scripts_dir,"pod/wavelet/params.r"))

##################################
### natty
##################################

### only wrf trials
trial_files <- list.files(wrf_les_cont_rel_trial_natty_dir, full.names = T)

### test
if(test){trial_files <- trial_files[1:2]}

cont_rel <- wave_detection(
	trial_files
	, nmin = nmin
	, detection_dir = wave_detection_dir
	, detection_file = 
		ifelse(
			test
			, "nat_detection_test.csv"
			, "nat_detection.csv" 
		)#edn ifelse
	, save_cont_rel = T
)#dedn wave detection

##################################
### enhanced
##################################

trial_files <- list.files(wrf_les_cont_rel_trial_enhanced_rho_dir, full.names = T)

### sample
trial_files <- sample(trial_files, wave_samples)

### test
if(test){trial_files <- list.files(wrf_les_cont_rel_trial_enhanced_test_dir, full.names = T)[42:125]}

cont_rel <- wave_detection(
	trial_files
	, nmin = nmin
	, detection_dir = wave_detection_dir
	, detection_file = 
		ifelse(
			test
			, "enhanced_detection_test.csv"
			, "enhanced_rho_detection.csv" 
		)#edn ifelse
	, save_cont_rel = T
)#dedn wave detection
