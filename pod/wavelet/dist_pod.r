
#########################
### examine PoD model of wave
#########################

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
source(file.path(scripts_dir,"pod/wavelet/params.r"))

list.files(wave_detection_dir)

###############################
### natty
###############################
nat_wave_detect_file <- file.path(wave_detection_dir, "nat_detection.csv")
nat_wave_cr <- read.csv(nat_wave_detect_file)

###############################
### enhanced
###############################

if(test){enh_wave_detect_file <- file.path(wave_detection_dir, "enhanced_detection_test.csv")}


enh_wave_detect_file <- file.path(wave_detection_dir, "enhanced_detection.csv")
# enh_rho_wave_detect_file <- file.path(wave_detection_dir, "enhanced_rho_detection.csv")


# enh_rho_wave_cr <- read.csv(enh_rho_wave_detect_file)
enh_wave_cr <- read.csv(enh_wave_detect_file)
plot_dd_hexmap(enh_wave_cr %>% filter(emiss_rate_kghr<2000),"wspd","emiss_rate_kghr")
# plot_dd_hexmap(enh_rho_wave_cr%>% filter(emiss_rate_kghr<2000) ,"wspd","emiss_rate_kghr")

# any(enh_wave_cr$DD & any(enh_wave_cr$emiss_rate_kghr==0))

# plot_dd_hexmap(enh_wave_cr %>% filter(emiss_rate_kghr<1000),"wspd","emiss_rate_kghr", nn = F) %>% save_png(fig_out_pod_wave_dir, "wave_heatmap.png")	


low <- enh_wave_cr %>% filter(emiss_rate_kghr>50,emiss_rate_kghr<75,gas_conc_noise<40,gas_conc_noise>20)
high <- enh_wave_cr %>% filter(emiss_rate_kghr>2000, area_pixel<2000)
xx <- low
xx <- rbind(low,high) 
xx <- xx[sample(1:dim(xx)[1],100),]

sum(xx$DD)/dim(xx)[1]
sum(xx$detached)/dim(xx)[1]

# xx <- xx[which(xx$emiss_rate_kghr == max(xx$emiss_rate_kghr)),]

lapply( list.files(pod_funcs_dir, full.names = T), source)
i <- 1 
for(i in 1:100){
	trial <- readRDS(xx[i,"trial_file"])
	scene <- reassemble_scene(trial)
	plot(scene)
	png(file.path( fig_out_pod_dir, paste0(i,".png")))
	plot(scene)
	dev.off()
}#edn for

	ww <- flip(t(wave(scene, hotspot =F)),"h")
	bb <- boundaries(ww)
	lines(as.polygons(bb))

	plot(scene)
	plot(ww)
	plot(flip(t(ww), "h"))
	plot(flip(ww))

	hh <- wave(scene, hotspot =T, nmin_hotspot = 1)
	DD <- detect_masked(hh, check_center =F)
	plot(c(scene,ww,hh)
	     , main = paste("PLUME:", i
			    ,"EMISS:", round(xx[i,"emiss_rate_kghr"])
			    ,"NOISE:", round(xx[i,"gas_conc_noise"])
			    ,"DETECT:", DD$detect
			    # ,"DETACH",DD$detached
	     )
	)
	Sys.sleep(2)
}#edn for

#### how often does wavelet find a clump in gauss noise? 
detects <- vector() 
detacheds <- vector() 
i <- 99
for(i in 1:100){
	print(i)
	values(scene) <- rnorm(length(values(scene)), mean = 0,sd = 30) +1800
	ww <- wave(scene)
	# plot(c(scene,ww))

	detect <- detect_masked(ww, check_center = T)
	detects <- c(detects,detect$detect)
	detacheds <- c(detacheds,detect$detached)
}##end fo

sum(detects)/length(detects)
sum(detacheds)/length(detacheds)

