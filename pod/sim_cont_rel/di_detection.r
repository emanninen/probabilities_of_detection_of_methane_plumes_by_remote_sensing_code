
#########################
### run di detection
### take controlled release scenes, run the di detection method 
### and get controlled release data with detect/no detect
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
### params
source(file.path(scripts_dir,"pod/sim_cont_rel/params.r"))

### generate trial scenes
# source(file.path(scripts_dir,"pod/sim_cont_rel/cont_rel.r"))

##########################
### natty run di detect
##########################
trial_files_natty <- list.files(wrf_les_cont_rel_trial_natty_dir, full.names = T)

cont_rel_natty <- di_detection(
	trial_files_natty
	, detection_dir = wrf_les_detection_dir
	, detection_file = "natty_di_detection.csv"
	, save_cont_rel = T
	, nmin = nmin
)#end cont_rel

##########################
### enhanced rung di detect
##########################

# trial_files_enhanced <- list.files(wrf_les_cont_rel_trial_enhanced_dir, full.names = T)
trial_files_enhanced <- list.files(wrf_les_cont_rel_trial_enhanced_rho_dir, full.names = T)

trial_files_enhanced <- sample(trial_files_enhanced, di_samples)

### test
if(test){trial_files_enhanced <- list.files(wrf_les_cont_rel_trial_enhanced_test_dir, full.names = T)}

cont_rel <- di_detection(
	trial_files_enhanced
	, nmin = nmin
	, hotspot = F
	, check_center = T
	, detection_dir = wrf_les_detection_dir
	, detection_file = 
		ifelse(
			test
			, "enhanced_di_detection_test.csv"
			# , "enhanced_di_detection.csv" 
			, "enhanced_rho_di_detection.csv"
		)#edn ifelse
	, save_cont_rel = T
)#dedn di detection

plot_dd_hexmap(cont_rel %>% filter(emiss_rate_kghr<1000),"wspd","emiss_rate_kghr")

if(F){
for(ii in 1:dim(cont_rel)[1]){
	print(ii)
	trial <- readRDS(cont_rel[ii,"trial_file"])
	scene <- reassemble_scene(trial)
	isol <- plume_isolate(scene, nmin=40, hotspot =F)
	isol_hot <- plume_isolate(scene, nmin=40)
	plot(c(scene,isol,isol_hot))
	Sys.sleep(2)
}#edn for(ii in 1:dim(cont_rel)[1]){
}#end iff 
