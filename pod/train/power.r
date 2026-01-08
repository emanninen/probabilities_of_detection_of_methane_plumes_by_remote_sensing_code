
#############################
### train power pod models like bridger
#############################

###  power model 
pow <- expo_log_reg(cont_rel_enh_di$train)
### This is worse than log reg
pow$coefficients

### ://stackoverflow.com/questions/45362548/nonlinear-logistic-regression-package-in-r
