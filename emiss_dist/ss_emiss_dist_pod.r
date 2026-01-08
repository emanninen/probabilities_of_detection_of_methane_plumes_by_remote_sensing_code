
###########################################
### calculate a analytical solution to the steady state distribution of point sources
### for different msat pod scenarios in Permian
###########################################

### dirs
while(!"dirs.r" %in% list.files()){setwd("..")};source("dirs.r")
### depends
lapply( list.files(emiss_dist_funcs_dir, full.names = T), source)
lapply( list.files(pod_funcs_dir, full.names = T), source)
library(dplyr)

set.seed(42)

### hmmm 
# source(file.path(emiss_dist_scripts_dir,"merge/main.r"))
# source(file.path(emiss_dist_scripts_dir,"merge/merge.r"))

### read in merged plume count data
merge_ct <- read.csv(file.path(plume_counts_dir,"merge_ct_pod.csv"))

### calculate steady state point source count 
### wiht different gcn pod

### only want 1 category for permian mair, which has different gcn for each scene
merge_ct$gas_conc_noise[merge_ct$system == "mair"] <- 41

ss_ps <- do.call(rbind,lapply(unique(merge_ct$gas_conc_noise),function(gcn){
	print(gcn)
	ct_gcn <- filter(merge_ct, gas_conc_noise == gcn)
	if(dim(ct_gcn)[1] == 0){return(NULL)}
	### dont analyze
	if(!any(ct_gcn$nn>0)){return(NULL)}
	### 
	ss_ps <- ana_ss_emiss_dist(
		ct_gcn
		, unique(merge_ct$emiss_rate_kghr_bin) 
		, pod_resample = T
	)#end anassemissdits
	ss_ps$gcn <- gcn
	return(ss_ps)
}))#edn for gcn

########################
### percent of total emissions
########################

### point source steady state total emiss
ss_ps_total_emiss <- ss_ps %>% 
	group_by(gcn) %>%
	dplyr::summarize( 
		 ### mair naive
		ps_emiss_total_naive = sum(ss_emiss_naive, na.rm = T)
		, ps_emiss_total_naive_lb = sum(ss_emiss_naive_lb, na.rm = T)
		, ps_emiss_total_naive_ub = sum(ss_emiss_naive_ub, na.rm = T)
		 ### mair pod
		, ps_emiss_total_pod = sum(ss_emiss_pod, na.rm = T)
		, ps_emiss_total_pod_lb = sum(ss_emiss_pod_lb, na.rm = T)
		, ps_emiss_total_pod_ub = sum(ss_emiss_pod_ub, na.rm = T)
	) %>% as.data.frame()#end sumamrzie

ss_ps <- merge(ss_ps, ss_ps_total_emiss, by = "gcn", all.x = T)

### read in area emissions
area_emiss <- read.csv(file.path(area_emiss_ss_dir, "area_emiss.csv"))
### only permian area emissions
permian_area_emiss <- area_emiss[area_emiss$basin=="Permian","area_emiss_gim"]

### calcualate steady state point source fraction of total emission rate
ss_ps$ss_rat_naive <- ss_ps$ss_emiss_naive / 
	(permian_area_emiss + ss_ps$ps_emiss_total_naive)
ss_ps$ss_rat_naive_lb <- ss_ps$ss_emiss_naive_lb / 
	(permian_area_emiss + ss_ps$ps_emiss_total_naive_lb)
ss_ps$ss_rat_naive_ub <- ss_ps$ss_emiss_naive_ub / 
	(permian_area_emiss + ss_ps$ps_emiss_total_naive_ub)

### pod
ss_ps$ss_rat_pod <- ss_ps$ss_emiss_pod / 
	(permian_area_emiss + ss_ps$ps_emiss_total_pod)
	# (permian_area_emiss + ss_ps$ss_emiss_pod)
ss_ps$ss_rat_pod_lb <- ss_ps$ss_emiss_pod_lb / 
	(permian_area_emiss + ss_ps$ps_emiss_total_pod_lb)
	# (permian_area_emiss + ss_ps$ss_emiss_pod)
ss_ps$ss_rat_pod_ub <- ss_ps$ss_emiss_pod_ub / 
	(permian_area_emiss + ss_ps$ps_emiss_total_pod_ub)
	# (permian_area_emiss + ss_ps$ss_emiss_pod)

write.csv(ss_ps, file.path(ss_emiss_dir, "ss_ps_pod.csv"), row.names = F)
