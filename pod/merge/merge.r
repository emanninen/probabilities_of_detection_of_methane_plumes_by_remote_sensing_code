
#########################
### merge cont rels from different experiments
#########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(dplyr)

#########################
#########################
### merge DI
#########################
#########################

#########################
### merge MAIR
#########################

list.files(mair_detection_dir)

### load in cont rel datta
di_mair_nat_cr <- read.csv(file.path(mair_detection_dir, "natty_di_detection.csv"))
di_mair_enh_cr <- read.csv(file.path(mair_detection_dir, "enhanced_di_detection.csv"))

### cat vars
di_mair_nat_cr$det_method <- "di"
di_mair_enh_cr$det_method <- "di"
di_mair_nat_cr$enhance <- F
di_mair_enh_cr$enhance <- T

#########################
### merge WRFLES
#########################

### load in cont rel datta
di_wrf_les_enh_cr <- read.csv(file.path(wrf_les_detection_dir, "enhanced_di_detection.csv"))
di_wrf_les_nat_cr <- read.csv(file.path(wrf_les_detection_dir, "natty_di_detection.csv"))
### cat vars
di_wrf_les_nat_cr$det_method <- "di"
di_wrf_les_enh_cr$det_method <- "di"
di_wrf_les_nat_cr$enhance <- F
di_wrf_les_enh_cr$enhance <- T

#########################
#########################
### merge wave
#########################
#########################

#########################
### merge natty
#########################

### load in cont rel datta
wave_nat_cr <- read.csv(file.path(wave_detection_dir, "nat_detection.csv"))
# wave_enh_cr <- read.csv(file.path(wave_detection_dir, "enhanced_detection.csv"))
wave_enh_cr <- read.csv(file.path(wave_detection_dir, "enhanced_detection.csv"))
### cat vars
wave_nat_cr$det_method <- "wave"
wave_enh_cr$det_method <- "wave"
wave_nat_cr$enhance <- F
wave_enh_cr$enhance <- T

#########################
### MERGE
#########################

### merge into one cr data
merge_cr <- rbind(
	di_mair_nat_cr 
	, di_mair_enh_cr
	, di_wrf_les_nat_cr
	, di_wrf_les_enh_cr
	### TODO: awkward
	# , wave_nat_cr
	, wave_enh_cr %>% mutate(check_center = F)
)#dnd rbind

### write df to disk
write.csv(
	merge_cr
	, file.path(merge_detection_dir ,  "detection.csv")
	,  row.names = F
)#end rwitecsv
