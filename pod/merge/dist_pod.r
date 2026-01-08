
#########################
### examine PoD model of merge
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotcli)

########################
### enhanced
########################
cont_rel <- read.csv(file.path(merge_detection_dir, "detection.csv"))
cont_rel <- cont_rel %>% filter(rf != "wrfles")
### keep to well sampled space
# cont_rel <- filter(cont_rel, emiss_rate_kghr<500, gas_conc_noise<30,gas_conc_noise>20,sqrt(area_pixel)<50)

	
# cont_rel <- cont_rel %>% filter(rf == "wrfles")
# plot_dd_xy(cont_rel,"ss","emiss_rate_kghr", "DD")
### This indicates an issue in the contrel detection 

plot_dd_hexmap(cont_rel,"emiss_rate_kghr","ss", nn = T)
plot_dd_hexmap(cont_rel,"ss","emiss_rate_kghr", nn = F)

# train <- train_pod(cont_rel)

train <- read.csv(file.path(pod_model_dir,"train.rds"))
# cont_rel <- pred_pod(cont_rel, train[["betas"]], train[["log_reg"]])
plot_pod_xx(cont_rel %>% filter(ss<100))

lapply( list.files(pod_funcs_dir, full.names = T), source)
gg <- redimension_grid(seq(10,500),wspd = 1.5)
hist100(gg$ss)

cont_rel <- cont_rel %>% filter(ss>50) 
hist100(cont_rel$ss)
hist100(cont_rel$emiss_rate_kghr)
hist100(cont_rel$area_pixel)
hist100(cont_rel$gas_conc_noise)
hist100(cont_rel$wspd)
gg <- pred_pod(gg, train[["betas"]], train[["log_reg"]])
plot_pod_xx(gg, "emiss_rate_kghr")

cr_samp <- cont_rel[sample(1:dim(cont_rel)[1],10),]

########################
### nattty
########################

nat_cr <- read.csv(file.path(merge_detection_dir, "natty_di_detection.csv"))
nat_cr <- pred_pod(nat_cr, train[["betas"]], train[["log_reg"]])
plot_pod_xx(nat_cr)
