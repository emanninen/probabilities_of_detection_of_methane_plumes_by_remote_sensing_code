
###################################
### Take a look at sampled condis from mair maps
###################################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)

mair_file <- list.files(thresh_dir,full.names =T )[22]
mair <- read.csv(mair_file)
