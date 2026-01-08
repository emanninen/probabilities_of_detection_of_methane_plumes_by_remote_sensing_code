
#########################################
### train pod models
#########################################

### merge together the cont rel data
source(file.path(pod_scripts_dir,"merge/merge.r"))

### split train test data
source(file.path(pod_scripts_dir,"train/split_train_test.r"))

### final model name
final_model <- "log_lr"

###  train di models
di <- train_pod(cont_rel_sim_di$train)[[final_model]]
###  select final moddel and save
saveRDS(di,file.path(pod_model_dir, "di.rds"))

### train models wave
wave <- train_pod(cont_rel_sim_wave$train)[[final_model]]
###  select final moddel and save
saveRDS(wave,file.path(pod_model_dir, "wave.rds"))
