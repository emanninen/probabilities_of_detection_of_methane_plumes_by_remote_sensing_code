
#########################
### investigate adding auto corr noise to wrf les cont rel
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
### params
source(file.path(scripts_dir,"pod/sim_cont_rel/params.r"))

#####################
### look at each mair cont rel scene 
### sample the none plume bits
### measure correlation of noise
#####################


cont_rel_natty <- read.csv(file.path(mair_detection_dir, "natty_di_detection.csv"))

ars <- vector()
for(ii in 1:dim(cont_rel_natty[1])){
	
	trial <- readRDS(cont_rel_natty[ii,"trial_file"])
	scene <- reassemble_scene(trial)
	plume_mask <- plume_isolate(scene,nmin = 100)
	scene_mask <- mask(scene,plume_mask, inverse = T)
	### keep even if no detect for realism
	scene_ar <- autocor(scene)[[1]]
	scene_ar <- autocor(scene_mask)[[1]]
	plot(scene, main = paste(scene_ar,cont_rel_natty[ii,"DD"]))
	Sys.sleep(1)
	ars <- c(ars,scene_ar )
}#edn forms
hist(ars,breaks = 50 )


####################
### compare different ways of making AR noise
### SCW wants me to do it w gauss filter: uhh
####################


if(F){
noise <- matrix(rnorm(
		dim(xch4_agg)[1] * dim(xch4_agg)[2]
		, 0
		, nn
	)#edn rnomr
	, dim(xch4_agg)[1] 
	, dim(xch4_agg)[2]
)#end atrxi

		noise <- neuRosim::spatialnoise(
				dim(xch4_agg)
				, nscan = 1
				, method = "corr"
				, sigma = nn 
				, rho = rho
			)[,,1] # end neuro sim

rho <- .5

plot(c(rast(xx),rast(noise)))
plot(rast(xx),rast(noise))
plot(rast(xch4_agg))
plot(rast(xch4_noise))
plot(rast(xch4_noise))

}
