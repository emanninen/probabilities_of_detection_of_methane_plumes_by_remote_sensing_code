
#########################
### run sim cont rel wrf les
### create controlled release trials out of wrf les plumes
### and mutate addtl plumes by changing emiss rate, pixel area, gcn
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
### params
source(file.path(scripts_dir,"pod/sim_cont_rel/params.r"))

### jus wrf les output files
# scene_files <- list.files(wrf_les_ideal_outputs_dir, full.names = T)
scene_files <- list.files(wrf_les_ideal_outputs_dir, full.names = T)
### Ju did the 50 m/s and 500(!) m/s wspd! 
scene_files <- grep("50", scene_files, value = T, invert = T)

### TODO: currently 1ms is missing xch4_1
scene_files <- scene_files[2:length(scene_files)] 

###  test
if(test){scene_files <- scene_files[c(7,10,14)]}

##########################
### natty convert wrf les plumes
##########################
### clear trial dir
file.remove(list.files(wrf_les_cont_rel_trial_natty_dir, full.names = T))

trial <- trial_scene_wrfles(
	scene_files = scene_files
	, trial_dir = wrf_les_cont_rel_trial_natty_dir 
	, time_index = seq(1,max_tt,round(max_tt/time_samples))
	, save = T
)#end sim_cont_rel_wrf

##########################
### enhanced convert wrf les plumes
##########################

### clear trial dir
if(test){
	file.remove(list.files(wrf_les_cont_rel_trial_enhanced_test_dir, full.names = T))}else{
	# file.remove(list.files(wrf_les_cont_rel_trial_enhanced_dir, full.names = T))
}#dnd else


print(paste("enh cont rel samples :"
	    , n_qq_samples *n_nn_samples *length(time_samples)*length(apf_samples)*length(scene_files)
	    ))

trial <- trial_scene_wrfles(
	scene_files = scene_files
	, enhance = T
	, target_emiss = qq_samples 
	, target_noise = nn_samples
	, gcn_rho = .7 
	, target_area_pixel_fact = apf_samples 
	, time_index = time_samples 
	, trial_dir = 
		ifelse(
			test
			, wrf_les_cont_rel_trial_enhanced_test_dir
			, wrf_les_cont_rel_trial_enhanced_rho_dir
			# , wrf_les_cont_rel_trial_enhanced_dir
		)#end if else
	, save = T
)#end sim_cont_rel_wrf
