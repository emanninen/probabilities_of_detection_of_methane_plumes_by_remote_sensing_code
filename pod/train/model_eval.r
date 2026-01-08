
#########################################
### compare different ways of modelling pod
#########################################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply(list.files(pod_funcs_dir, full.names = T), source)
library(ggplot2)
### split test train
source(file.path(pod_scripts_dir,"train/split_train_test.r"))

##############################
### DI
##############################

# plot_dd_hexmap(cont_rel_sim_di$train %>% filter(emiss_rate_kghr<1000),"wspd","emiss_rate_kghr")
### train models
train_di <- train_pod(cont_rel_sim_di$train, snr_cov_scale = 5)
names(train_di)
train_di[["conrad"]] <- "conrad"

### predict on simulated test data w/ pred models
test_di <- pred_pod(cont_rel_sim_di$test,train_di)

### invesitage model tested on enhanced test data 
# plot_pod_xx(test_di, pod_varnames = names(train_di),x = "ss", x_binwidth = 4)
# test_di$log_ss <- log(test_di$ss)
# plot_pod_xx(test_di, pod_varnames = names(train_di),x = "log_ss", x_binwidth = 1e-3)
# plot_pod_xx(test_di, pod_varnames = names(train_di),x = "emiss_rate_kghr", x_binwidth = 100) + labs(x = "Emission Rate (kg/hr)")

### predict simulation trained pod model on test set of simulated scenes
# plot_roc(lapply_rocdf(test_di, pod_varnames = names(train_di)))

### predict simulation trained pod model on enhanced MAIR controlled release data
# pred_di <- pred_pod(cont_rel_di,train_di)
# plot_di_roc <- plot_roc(lapply_rocdf(pred_di, pod_varnames = names(train_di))) +
# 	labs(
# 		title = "Receiver Operator Curves for Divergence Integral Pd Models"
# 	) + # end labs



plot_di_roc <- plot_roc(
	lapply_rocdf(
	insightm
	, pod_varnames = c("lr","log_lr","conrad") 
	, length.out = 5000
	)#end lappl
)#dnd plot 
# +
# 	labs(
# 		title = "Receiver Operator Curves for Divergence Integral Pd Models"
# 	) # end labs
	

##############################
### WAVE
##############################

plot_dd_hexmap(cont_rel_sim_wave$train %>% filter(emiss_rate_kghr<1000),"wspd","emiss_rate_kghr")
### train models
train_wave <- train_pod(cont_rel_sim_wave$train, snr_cov_scale = 5)
names(train_wave)

### prewavect on enhanced test data w/ pred models
test_wave <- pred_pod(cont_rel_sim_wave$test,train_wave)
test_wave <- test_wave %>% filter(emiss_rate_kghr<1000)

### invesitage model tested on enhanced test data 
plot_pod_xx(test_wave, pod_varnames = names(train_wave),x = "ss", x_binwidth = 4)
test_wave$log_ss <- log(test_wave$ss)
plot_pod_xx(test_wave, pod_varnames = names(train_wave),x = "log_ss", x_binwidth = 1e-2)

plot_pod_xx(test_wave, pod_varnames = names(train_wave),x = "emiss_rate_kghr", x_binwidth = 100) +
	labs(x = "Emission Rate (kg/hr)")
plot_roc(lapply_rocdf(test_wave, pod_varnames = names(train_wave)))

### prewavect trained pod model on natty data
### TODO: somehow this isnt working in pod merge
cont_rel_nat_wave <- read.csv(file.path(wave_detection_dir, "nat_detection.csv"))
nat_wave <- pred_pod(cont_rel_nat_wave,train_wave)
