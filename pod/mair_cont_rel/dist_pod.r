
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

### hmmm
# source(file.path(scripts_dir,"pod/mair_cont_rel/run_di_detection.r"))

### read in controlled release results
cont_rel_natty <- read.csv(file.path(mair_detection_dir, "natty_di_detection.csv"))
plot_dd_hexmap(cont_rel_natty,"wspd","emiss_rate_kghr")

### read in controlled release results
cont_rel_enhanced <- read.csv( file.path(mair_detection_dir, ifelse(test,"enhanced_di_detection_test.csv","enhanced_di_detection.csv")))

unique(cont_rel_enhanced$wspd)
dim(cont_rel_enhanced)

plot_dd_hexmap(cont_rel_enhanced %>% filter(emiss_rate_kghr<2000,gas_conc_noise>30),"wspd","emiss_rate_kghr",nn =F,bins = 30)

