
#########################
### look at merging cont rel osurce
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(ggplot2)
library(cowplot)
library(dplyr)

list.files(mair_detection_dir)
mair_enh_cr <- read.csv(file.path(mair_detection_dir, "enhanced_di_detection.csv"))
wrf_les_enh_cr <- read.csv(file.path(wrf_les_detection_dir, "enhanced_di_detection.csv"))
names(wrf_les_enh_cr)%in% names(mair_enh_cr)
names(mair_enh_cr)%in% names(wrf_les_enh_cr)
names(mair_enh_cr)

merge_cr <- rbind(mair_enh_cr, wrf_les_enh_cr)

hist_var <- function(cr,x){
	plot_ss_hist <- ggplot(cr) + 
		geom_histogram(aes(x = .data[[x]],fill = rf)) + 
		theme_classic() +
		theme(text = element_text(size = 30))
	return(plot_ss_hist)
}#end func hist var

plot_grid(
	# hist_var(merge_cr %>% filter(rf != "wrfles"), "ss")
	hist_var(merge_cr, "ss")
	, hist_var(merge_cr, "emiss_rate_kghr")
	, hist_var(merge_cr, "emiss_rate_ppbm2s")
	# , hist_var(merge_cr, "p_surf")
	# , hist_var(merge_cr, "area_pixel")
	# , hist_var(merge_cr, "gas_conc_noise")
	# , hist_var(merge_cr, "wspd")
)#dend plot gricd

mair <- readRDS(list.files(mair_cont_rel_trial_enhanced_dir,full.names = T)[1])
sim <- readRDS(list.files(wrf_les_cont_rel_trial_enhanced_dir,full.names = T)[1])
mair$p_surf
sim$p_surf

