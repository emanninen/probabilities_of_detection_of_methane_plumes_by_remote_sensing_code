
#############################
### examine the PoD distributions produced by controlled release data
### sand boxy
#############################

### reproducibility
rm()
set.seed(42)
### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotcli)
library(terra)
### params
source(file.path(scripts_dir,"pod/sim_cont_rel/params.r"))

### hmmm
# source(file.path(scripts_dir,"pod/sim_cont_rel/di_detection.r"))

if(F){
#############################
### NATTY 
#############################
### read in controlled release results
cont_rel_natty <- read.csv(file.path(wrf_les_detection_dir, "natty_di_detection.csv"))
### set up p s grid over which to plot dist
grid_natty <- p_s_grid(min_ss =0, max_ss = .9*max(cont_rel_natty$ss))
### calculate beta binom post dist from cont rel data
stats_natty <- stats_fpod(cont_rel_natty, grid_natty, snr_cov_scale = 5e-3)
### quick plots
pp <- ggplot(cont_rel_natty,aes(x = ss, y = wspd, color = DD)) + geom_point()
ggplotcli(pp, braille = F, plot_width = 20)
}#end if f

#############################
### ENHANCED
#############################
### read in controlled release results
enh_detection_file <- ifelse(
	test
	, file.path(wrf_les_detection_dir, "enhanced_di_detection_test.csv")
	# , file.path(wrf_les_detection_dir, "enhanced_di_detection.csv")
	, file.path(wrf_les_detection_dir, "enhanced_rho_di_detection.csv")
)#edn ifsle
cont_rel_enhanced <- read.csv(enh_detection_file)
plot_dd_hexmap(cont_rel_enhanced,"wspd","emiss_rate_kghr") 

plot_dd_hexmap(cont_rel_enhanced%>% filter(emiss_rate_kghr<1000,gas_conc_noise<35),"wspd","emiss_rate_kghr") 
### none of the zero emiss rate were detected!
any(cont_rel_enhanced[cont_rel_enhanced$emiss_rate_kghr == 0,"DD"])

### quick plots
# pp <- ggplot(cont_rel_enhanced,aes(x = ss, y = wspd, color = DD)) + geom_point()
# ggplotcli(pp, braille = F, plot_width = 20, plot_height = 20)
cont_rel_enhanced$log_ss <- log(cont_rel_enhanced$ss)
plot_dd_line( cont_rel_enhanced %>% filter(emiss_rate_kghr<1000))#end plotdd line
plot_dd_line( cont_rel_enhanced %>% filter(emiss_rate_kghr<1000),"log_ss")#end plotdd line

plot_dd_hexmap(cont_rel_enhanced%>% filter(emiss_rate_kghr<1000,gas_conc_noise<35),"wspd","emiss_rate_kghr") 
%>% save_png(fig_out_pod_di_dir, "di_heatmap.png")	

