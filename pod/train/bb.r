
###################################
### train beta binomial models of pod
###################################

### calculate beta binom model and log reg pod model on trainign data
### TODO: split out models if you add more models
train <- train_pod(train_cr, snr_cov_scale)
