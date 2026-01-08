
#########################
### params for sim cont rel
#########################

###reproduce
rm()
set.seed(42)

### test?
test <- T

#########################
### cont_rel.r
#########################

### this kind of lines up with the space fo the cont rel exp
min_qq <- 1
max_qq <- 150
mean_qq <- 6
sd_qq <- 1.1
n_qq_samples <- 50
if(test){n_qq_samples <- 5}
qq_samples <- rlnorm(n_qq_samples,mean_qq,sd_qq)
which(qq_samples<100)

min_nn <- 10
max_nn <- 60
mu_nn <- 26 #ppb
sd_nn <- 6 #ppb
n_nn_samples <- 30 
if(test){n_nn_samples <- 5}
nn_samples <-  rnorm(n_nn_samples,mu_nn,sd_nn)

min_ap_fact <- 1
max_ap_fact<- 2
apf_samples <- min_ap_fact:max_ap_fact

# n_samples_per_cr <- 200
time_samples <- 10
max_tt <- 500 


time_samples <- seq(1,max_tt,round(max_tt/time_samples))
if(test){time_samples <- 5}

### test
# if(test){n_samples_per_cr <- 100}


#########################
### di_detection
#########################

nmin <- 40
di_samples <- 50e3
