
########################
### merge controlled release results
### copy and pasted from Ju's paper
### to a single file
########################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
library(lubridate)

list.files(mair_cont_rel_meter_emiss_rate_dir)
mair_blind <- read.table(file.path(mair_cont_rel_meter_emiss_rate_dir, "S2_blind_mair.txt"), header = T)
mair_meter <- read.table(file.path(mair_cont_rel_meter_emiss_rate_dir, "S5_mair_meter.txt"), header = T)
maire <- read.table(file.path(mair_cont_rel_meter_emiss_rate_dir, "S6_maire.txt"), header = T)

meter <- read.csv(file.path(mair_cont_rel_meter_emiss_rate_dir, "meterDF_MAIR_unblindedToMAIR.csv"))
matched_meter <- read.csv(file.path(mair_cont_rel_meter_emiss_rate_dir, "matchedDF_MAIR_unblindedToMAIR.csv"))

mair <- mair_meter
### bruh... I fixed the date ...
mair_meter$cr_kgh_CH4_mean90 == mair_blind$cr_kgh_CH4_mean90
mair$blinded <- mair_blind$blinded
mair$utc <- dmy_hms(mair$utc)

maire[,c("mair_l", "mair_u")] <- NULL
colnames(maire)[colnames(maire) == "mair"] <- "blinded"
maire$utc <- dmy_hms(maire$utc)
maire$Winds <- NA 
mair

meter <- rbind(mair, maire)

write.csv(meter, file.path(mair_cont_rel_meter_emiss_rate_dir, "ju_meter.csv"))
