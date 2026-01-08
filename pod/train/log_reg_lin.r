
### TODO: read in cont rel data

### train linear pod model
lr_lin <- glm(
	DD~emiss_rate_kghr+area_pixel+gas_conc_noise+wspd
	, data = train_cr
	, family = "binomial"
)#end glm

### TODO: save model
