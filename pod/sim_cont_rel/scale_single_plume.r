
###########################
### take a single wrf lesplume 
### scale emission rate, and noise, see detection
###########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(ggplot2)
library(dplyr)

### read in controlled release results
enh_detection_file <-  file.path(wrf_les_detection_dir, "enhanced_di_detection_test.csv")
cont_rel_enhanced <- read.csv(enh_detection_file)

plot_dd_hexmap(cont_rel_enhanced %>% filter(ss<5),'wspd', "emiss_rate_kghr")

sf <- unique(cont_rel_enhanced$scene_file)[2]

unique(cont_rel_enhanced$tt)

tti <- 57
cr <- cont_rel_enhanced %>%
	filter(
		# scene_file == sf
		, tt == tti
	        , !DD
	)#end filter

unique(cr$scene_file)
unique(cr$tt)
ii <- 1

for(ii in 1:dim(cr)[1]){
	cc <- cr[ii,]
	tf <- readRDS(cc$trial_file)
	ss <- reassemble_scene(tf)
	plot(ss)
	plot( plume_isolate(ss,nmin = 40))
	Sys.sleep(.5)
}#end

plot_dd_hexmap(cr,'', "emiss_rate_kghr")
cr$DD
cr <- cr[200,]
trial <- readRDS(cr$trial_file)
scene <- reassemble_scene(trial)
xx <- plume_isolate(scene, nmin =100)
