
#########################
### params for sim cont rel
#########################

set.seed(42)

test <- F

#########################
### cont_rel.r
#########################

rf_dirs <- list(
	mair_cont_rel_rf04
	, mair_cont_rel_rf05
	, mair_cont_rel_rf01e
	, mair_cont_rel_rf03e
)#end list

rfs <- c(
	"rf01e"
	, "rf03e"
	, "rf04"
	, "rf05"
)#edn list

cont_rel_locations <- list( 
	rf04 = c(-102.301, 32.053)
	, rf05 = c(-102.301, 32.053)
	, rf01e = c(-111.785, 32.821)
	, rf03e = c(-111.785, 32.821)
)#end list

### this kind of lines up with the space fo the cont rel exp
mean_qq <- 6
sd_qq <- 1.1
n_qq_samples <- 120
if(test){n_qq_samples <- 5}
qq_samples <- rlnorm(n_qq_samples,mean_qq,sd_qq)

mu_nn <- 26 #ppb
sd_nn <- 6 #ppb
n_nn_samples <- 20 
if(test){n_nn_samples <- 5}
nn_samples <-  rnorm(n_nn_samples,mu_nn,sd_nn)

min_ap_fact <- 1
max_ap_fact<- 2
apf_samples <- min_ap_fact:max_ap_fact
### on a odyssyey interactive node, takes ~.1 s per sample
### ~100 natty mair cont rel scenes
### time for simulated cont rel ~10s * n samples

#########################
### di_detection
#########################

nmin <- 100
di_samples <- 50e3
