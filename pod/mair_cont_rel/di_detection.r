
##########################
### run only the plume finding section of maryanns DI algo to look for a plume
### in rfo4, rf05 controlled release trial scene
##########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
### params
source(file.path(scripts_dir,"pod/mair_cont_rel/params.r"))

### hmmm
source(file.path(scripts_dir,"pod/mair_cont_rel/cont_rel.r"))

#########################
### natty
#########################

trial_files <- list.files(mair_cont_rel_trial_natty_dir, full.names = T)
cont_rel <- di_detection(trial_files
			 , detection_dir = mair_detection_dir
			 , detection_file = "natty_di_detection.csv"
			 , save_cont_rel = T
)#end di detect

#########################
### enhanced
#########################

trial_files <- list.files(mair_cont_rel_trial_enhanced_dir, full.names = T)
if(test){trial_files <- list.files(mair_cont_rel_trial_enhanced_test_dir, full.names = T)}

trial_files <- sample(trial_files, di_samples)

cont_rel <- di_detection(
	trial_files
	, detection_dir = mair_detection_dir
	, detection_file = 
		ifelse(
			test
			, "enhanced_di_detection_test.csv"
			, "enhanced_di_detection.csv" 
		)#edn ifelse
	, save_cont_rel = T
)#end di_detection
