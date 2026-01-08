
##################################
### split test and train data sets off the merge
##################################

### reproduce
rm()
set.seed(42)
### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply(list.files(pod_funcs_dir, full.names = T), source)
library(dplyr)

### read in merged cont rel data
cont_rel <- read.csv(file.path(merge_detection_dir, "detection.csv"))

#########
### sim for traingin
#########
cont_rel_sim <- cont_rel %>% filter(rf == "wrfles")
# cont_rel <- cont_rel %>% filter(rf != "wrfles")

### dont bother splitting out enhanced vs natty for simulation training
###  split out wave/di
cont_rel_sim_di <- cont_rel_sim %>% filter(det_method == "di")
cont_rel_sim_wave <- cont_rel_sim %>% filter(det_method == "wave")

cont_rel_sim_di <- split_test_train(cont_rel_sim_di)
cont_rel_sim_wave <- split_test_train(cont_rel_sim_wave)

#########
### cont rela
#########

cont_rel <- cont_rel %>% filter(rf != "wrfles")

### dont bother splitting out enhanced vs natty for simulation training
###  split out wave/di
cont_rel_di <- cont_rel %>% filter(det_method == "di")
cont_rel_wave <- cont_rel %>% filter(det_method == "wave")

### no test train split for cont rel b/c only used for validation
