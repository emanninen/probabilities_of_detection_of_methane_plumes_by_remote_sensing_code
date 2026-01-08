
################################
### plot example wrf les plumes
################################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
plumes_dir <- "/storage/research/emiss_dist/pod/cont_rel/sim/trials/enhanced_rho"
### depends
lapply(list.files(file.path(phd_dir,"r_utility"), full.names = T), source)
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(stringr)
library(ggplot2)
library(terra)
library(viridis)

set.seed(41)

plume_files <- list.files(plumes_dir, full.names = T)

### examine 1 to see components in filenmae
xx <- readRDS(plume_files[1]) 
names(xx)
xx$prior_wspd
xx$emiss_rate_kghr

### get filenames
plume_names <- str_split_i(plume_files, "/", -1)

### get components
plumes <- data.frame(
	plume_files = plume_files
	, plume_names = plume_names
	, area_pixel = as.numeric(str_split_i(plume_names, "_", -2))
	, gcn = as.numeric(str_split_i(plume_names, "_", -3))
	, wspd = as.numeric(str_split_i(plume_names, "_", -4))
	, qq = as.numeric(str_split_i(plume_names, "_", -5))
)#edn df

### select trials
### apparently plotting different resolution ones is too hard
### so 
plumes <- plumes[plumes$area_pixel==1089,]
plumes_sample <- plumes[sample(1:dim(plumes),9),]
### read in trials
trials <- lapply(plumes_sample$plume_files, readRDS)
scenes <- lapply(trials, reassemble_scene)

standard_scene <- scenes[[1]] 
standard_scene_df <- as.data.frame(scenes[[1]] , xy = T)

ss <-4 #test
scenes_df <- data.frame() 
for(ss in 1:length(scenes)){
	print(ss)
	scene <- scenes[[ss]]
	ext(scene)[2] <51
	if(ext(scene)[2]<51){
		ext(scene) <- ext(standard_scene)
	}#end if
	scene <- as.data.frame(scene,xy = T)
	scene$tag <- as.character(ss)
	scenes_df <- rbind(scenes_df, scene)
}#end forssscenes_df

### facet labels
facet_labels <- paste(
	      "Emiss. Rate:"
	      , round(plumes_sample$qq)
	      , "|"
	      , "Noise:"
	      , round(plumes_sample$gcn)
	      , "|"
	      , "WSPD:"
	      , round(plumes_sample$wspd,1)
)#enc pastee
names(facet_labels) <- 1:9
# facet_labels <- paste("Emiss. Rate:", plumes_sample$qq, "Noise:", plumes_sample$gcn, "WSPD:", plumes_sample$wspd)
### plot trial plumes on a grid
plot_trials <- ggplot(scenes_df) +
	geom_tile(aes(x = x, y = y, fill = lyr.1)) +
	# facet_wrap(~tag, scales = "free", labeller = as_labeller(facet_labels)) +
	facet_wrap(~tag, labeller = labeller(tag = facet_labels)) +
	labs(
		title =  "Example Simulated Plumes"
		, fill = expression(paste(XCH[4]),"(ppm)")

	) +
	scale_fill_viridis(option = "turbo", end = .9) +
	ethans_theme() +
	theme(strip.text = element_text(size = 25)) +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
	theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
	theme(legend.position = "bottom") +
	theme(legend.key.width = unit(5,"cm"))

print(plot_trials)

save_png(plot_trials, fig_out_pod_dir, "trial_scenes.png", width = 1600)
