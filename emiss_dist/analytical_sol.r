
#######################################
### analytical soln function
#######################################

########################
### ana_ss_emiss_dist 
### calculate analytical steady state emiss distribtuion population
### for a single basin/campaign
########################
ana_ss_emiss_dist <- function(
	ct
	, emiss_rate_kghr_bin_full
	, n_boot = 1e3
	, pod_resample = F
){

	basin_area <- unique(ct$basin_area)

	ss_ps <- data.frame(
		emiss_rate_kghr_bin = unique(ct$emiss_rate_kghr_bin)[order(unique(ct$emiss_rate_kghr_bin))]
	)#end dataframe

	### area scaled by pod is bottom of ss fractionform 
	campaign_area <- ct %>% 
		group_by(emiss_rate_kghr_bin) %>%
		dplyr::summarize(
			total_area_km2 = sum(mair_area_km2)
			, total_area_pod_km2 = sum(mair_area_km2*pod)
		) %>% as.data.frame()
	### order by emiss
	campaign_area <- campaign_area[order(campaign_area$emiss_rate_kghr_bin),]

	ii <- 1000 # test
	#################
	### bootstrap 
	#################
	ps <- sapply(1:n_boot, function(ii){
		### boot sample
		ct_boot <- ct[sample(1:dim(ct)[1], replace = T),]
		### did msat see each one that mair saw?
		if(pod_resample){
				sampled_nn <- as.vector(apply(
					ct_boot[,c("nn","pod")]
					, 1
					, function(cc){rbinom(1,cc[["nn"]],cc[["pod"]])}
				))#end apply
				ct_boot$nn <- sampled_nn
		}#dend if pod resamepl
		### reform to be counts for each emiss bin
		ps_boot <- nn_bin(ct_boot, emiss_rate_kghr_bin_full)
		return(ps_boot$ps)

	}) #end sapplye

	######################
	### Summarize boots
	### scale point source seen to steady state MLE point sources
	######################

	### Boot CI's naive
	ss_ps$ss_ps_naive <- apply(ps, 1,mean, na.rm = T) *
		(basin_area/campaign_area$total_area_km2)
	ss_ps$ss_ps_naive_lb <- apply(ps, 1, quantile, .025, na.rm = T) *
		(basin_area/campaign_area$total_area_km2)
	ss_ps$ss_ps_naive_ub <- apply(ps, 1, quantile, .975, na.rm = T) *
		(basin_area/campaign_area$total_area_km2)

	### Boot CI's pod
	ss_ps$ss_ps_pod <- apply(ps, 1, mean, na.rm = T) *
		(basin_area/campaign_area$total_area_pod_km2)
	ss_ps$ss_ps_pod_lb <- apply(ps, 1, quantile, .025, na.rm = T) *
		(basin_area/campaign_area$total_area_pod_km2)
	ss_ps$ss_ps_pod_ub <- apply(ps, 1, quantile, .975, na.rm = T) *
		(basin_area/campaign_area$total_area_pod_km2)

	######################
	### calcualte emiss rate
	######################

	### calculate steady state emiss rate naive
	ss_ps$ss_emiss_naive <- ss_ps$ss_ps_naive*ss_ps$emiss_rate_kghr_bin
	ss_ps$ss_emiss_naive_lb <- ss_ps$ss_ps_naive_lb*ss_ps$emiss_rate_kghr_bin
	ss_ps$ss_emiss_naive_ub <- ss_ps$ss_ps_naive_ub*ss_ps$emiss_rate_kghr_bin

	### calculate steady state emiss rate pod
	ss_ps$ss_emiss_pod <- ss_ps$ss_ps_pod*ss_ps$emiss_rate_kghr_bin
	ss_ps$ss_emiss_pod_lb <- ss_ps$ss_ps_pod_lb*ss_ps$emiss_rate_kghr_bin
	ss_ps$ss_emiss_pod_ub <- ss_ps$ss_ps_pod_ub*ss_ps$emiss_rate_kghr_bin

	return(ss_ps)
}#end fucn

#############
### ss_emiss_bin 
### steady state ps in each emiss rate bbin	
#############
nn_bin <- function(ct, emiss_rate_kghr_bin_full){
	ss_ps <- ct %>% 
		group_by(emiss_rate_kghr_bin) %>%
		dplyr::summarize(
			ps = sum(nn)
		) %>% as.data.frame()

	### if an emiss rate bin is not represented, add it w/ zero ps
	if(any(!emiss_rate_kghr_bin_full %in% ct$emiss_rate_kghr_bin)){
		missing_emiss_bins <- unique(emiss_rate_kghr_bin_full[
			which(
			      !emiss_rate_kghr_bin_full %in% 
			      ct$emiss_rate_kghr_bin
			)#end whcih
		])#ddn [)

		missing_emiss_bins <- cbind(missing_emiss_bins,0)
					    ### guh
					    # rep(0,
						# dim(ss_ps)[2]-1
					    # )
		# )
		colnames(missing_emiss_bins) <- c(
			"emiss_rate_kghr_bin"
			,"ps"
		)#ednd c 
		ss_ps <- rbind(ss_ps, missing_emiss_bins)
	}#edn if
	return(ss_ps)
}#edn func
