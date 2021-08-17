

# Sampled trajectories:
if(run_pars$trapfit==TRUE){
	fig_ct_fit <- plot_ct_fit_trap(params_df, global_pars, indiv_data, ctalpha=0.01, xlim=c(-20,30), ntraces=100)
	} else{
	fig_ct_fit <- plot_ct_fit(params_df, global_pars, indiv_data, ctalpha=0.01, xlim=c(-20,30), ntraces=100)		
	}

if(run_pars$trapfit==TRUE){
	fig_ct_fit_green <- plot_ct_fit_trap(params_df, global_pars, indiv_data, ctalpha=0.01, specialcol="green", xlim=c(-20,30), ntraces=100)
	} else{
	fig_ct_fit_green <- plot_ct_fit(params_df, global_pars, indiv_data, ctalpha=0.01, specialcol="green", xlim=c(-20,30), ntraces=100)		
	}

if(run_pars$trapfit==TRUE){
	fig_ct_fit_black <- plot_ct_fit_trap(params_df, global_pars, indiv_data, ctalpha=0.01, specialcol="black", xlim=c(-20,30), ntraces=100)
	} else{
	fig_ct_fit_black <- plot_ct_fit(params_df, global_pars, indiv_data, ctalpha=0.01, specialcol="black", xlim=c(-20,30), ntraces=100)		
	}

# # Sampled trajectories:
# fig_ct_fit_special_fulldata <- plot_ct_fit_special_fulldata(params_df, global_pars, ct_dat_refined, ctalpha=0.01, ntraces=100)

# # Sampled trajectories:
# fig_ct_fit_special_noid <- plot_ct_fit_special_noid(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=100)

# # Sampled trajectories:
# fig_ct_fit_special_fulldata_noid <- plot_ct_fit_special_fulldata_noid(params_df, global_pars, ct_dat_refined, ctalpha=0.01, ntraces=100)


if(run_pars$trapfit==TRUE){
	meanvalsindiv <- params_df %>% 
		mutate(infdur=wp+wx+wr) %>%
		group_by(id) %>% 
		summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wx=mean(wx), wr=mean(wr), infdur=mean(infdur), special=first(special))
	} else{
	meanvalsindiv <- params_df %>% 
		mutate(infdur=wp+wr) %>%
		group_by(id) %>% 
		summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr), infdur=mean(infdur), special=first(special))
	}


jitterfactor <- 0.01
fig_dpmean <- shared_params_df %>% 
	select(dpmeanB_trans, dpmeanW_trans) %>% 
	rename(dpmeanB=dpmeanB_trans, dpmeanW=dpmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=global_pars[["lod"]]-value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(10,NA)) + 
		scale_color_manual(values=c("dpmeanB"="red","dpmeanW"="blue","1"="red","0"="blue"), labels=c("dpmeanB"="Variant", "dpmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("dpmeanB"="red","dpmeanW"="blue"), labels=c("dpmeanB"="Variant", "dpmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		# geom_point(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=jitterfactor*special, col=factor(special)))  + 
		# geom_vline(data=meanvalsindiv, aes(xintercept=global_pars[["lod"]]-dp, col=factor(special)), alpha=0.2)  + 
		labs(x="Mean peak Ct", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_dpmean_withprior <- shared_params_df %>% 
	select(dpmeanB_trans, dpmeanW_trans) %>% 
	rename(dpmeanB=dpmeanB_trans, dpmeanW=dpmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=global_pars[["lod"]]-value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(10,NA)) + 
		scale_color_manual(values=c("dpmeanB"="red","dpmeanW"="blue","1"="red","0"="blue"), labels=c("dpmeanB"="Variant", "dpmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("dpmeanB"="red","dpmeanW"="blue"), labels=c("dpmeanB"="Variant", "dpmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		# geom_point(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=jitterfactor*special, col=factor(special)))  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$dpmean_prior, sd=prior_pars$dpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		# geom_vline(data=meanvalsindiv, aes(xintercept=global_pars[["lod"]]-dp, col=factor(special)), alpha=0.2)  + 
		labs(x="Mean peak Ct", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_gemlmean <- shared_params_df %>% 
	select(dpmeanB_trans, dpmeanW_trans) %>% 
	rename(dpmeanB=dpmeanB_trans, dpmeanW=dpmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	mutate(value=10^convert_Ct_logGEML(global_pars[["lod"]]-value)) %>%
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)), breaks=c(10^6, 10^7, 10^8, 10^9, 10^10)) + 
		# scale_x_continuous(limits=c(10,NA)) + 
		scale_color_manual(values=c("dpmeanB"="red","dpmeanW"="blue","1"="red","0"="blue"), labels=c("dpmeanB"="Variant", "dpmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("dpmeanB"="red","dpmeanW"="blue"), labels=c("dpmeanB"="Variant", "dpmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=10^convert_Ct_logGEML(global_pars[["lod"]]-dp), y=jitterfactor*special*5, col=factor(special)), alpha=0.4, width=0, height=0.003*5)  + 
		# geom_point(data=meanvalsindiv, aes(x=10^convert_Ct_logGEML(global_pars[["lod"]]-dp), y=jitterfactor*special*5, col=factor(special)))  + 
		# geom_vline(data=meanvalsindiv, aes(xintercept=global_pars[["lod"]]-dp, col=factor(special)), alpha=0.2)  + 
		labs(x="Mean peak RNA copies per ml", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off


fig_wpmean <- shared_params_df %>% 
	select(wpmeanB_trans, wpmeanW_trans) %>% 
	rename(wpmeanB=wpmeanB_trans, wpmeanW=wpmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		# geom_histogram(aes(col=name, fill=name, y=..density..), alpha=0.2, bins=100, position="identity") + 
		scale_x_continuous(limits=c(0,NA)) + 
		# scale_color_manual(values=c("wpmeanB"="red","wpmeanW"="blue"), labels=c("wpmeanB"="Variant", "wpmeanW"="Non-variant")) + 
		scale_color_manual(values=c("wpmeanB"="red","wpmeanW"="blue","1"="red","0"="blue"), labels=c("wpmeanB"="Variant", "wpmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("wpmeanB"="red","wpmeanW"="blue"), labels=c("wpmeanB"="Variant", "wpmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=wp, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		# geom_point(data=meanvalsindiv, aes(x=wp, y=jitterfactor*special, col=factor(special)))  + 
		labs(x="Mean proliferation stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_wpmean_withprior <- shared_params_df %>% 
	# select(wpmeanB, wpmeanW) %>% 
	select(wpmeanB_trans, wpmeanW_trans) %>% 
	rename(wpmeanB=wpmeanB_trans, wpmeanW=wpmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		# scale_color_manual(values=c("wpmeanB"="red","wpmeanW"="blue"), labels=c("wpmeanB"="Variant", "wpmeanW"="Non-variant")) + 
		scale_color_manual(values=c("wpmeanB"="red","wpmeanW"="blue","1"="red","0"="blue"), labels=c("wpmeanB"="Variant", "wpmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("wpmeanB"="red","wpmeanW"="blue"), labels=c("wpmeanB"="Variant", "wpmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=wp, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$wpmean_prior, sd=prior_pars$wpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		labs(x="Mean proliferation stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

if(run_pars$trapfit==TRUE){
fig_wxmean <- shared_params_df %>% 
	# select(wxmeanB, wxmeanW) %>% 
	select(wxmeanB_trans, wxmeanW_trans) %>% 
	rename(wxmeanB=wxmeanB_trans, wxmeanW=wxmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("wxmeanB"="red","wxmeanW"="blue","1"="red","0"="blue"), labels=c("wxmeanB"="Variant", "wxmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("wxmeanB"="red","wxmeanW"="blue"), labels=c("wxmeanB"="Variant", "wxmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=wx, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		labs(x="Mean plateau stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_wxmean_withprior <- shared_params_df %>% 
	# select(wxmeanB, wxmeanW) %>% 
	select(wxmeanB_trans, wxmeanW_trans) %>% 
	rename(wxmeanB=wxmeanB_trans, wxmeanW=wxmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("wxmeanB"="red","wxmeanW"="blue","1"="red","0"="blue"), labels=c("wxmeanB"="Variant", "wxmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("wxmeanB"="red","wxmeanW"="blue"), labels=c("wxmeanB"="Variant", "wxmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=wx, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$wxmean_prior, sd=prior_pars$wxsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		labs(x="Mean plateau stage duration (days)", y="Density") + 
		theme_minimal() + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off
}

fig_wrmean <- shared_params_df %>% 
	# select(wrmeanB, wrmeanW) %>% 
	select(wrmeanB_trans, wrmeanW_trans) %>% 
	rename(wrmeanB=wrmeanB_trans, wrmeanW=wrmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		# geom_histogram(aes(col=name, fill=name, y=..density..), alpha=0.2, bins=100, position="identity") + 
		scale_x_continuous(limits=c(0,NA)) + 
		# scale_color_manual(values=c("wrmeanB"="red","wrmeanW"="blue"), labels=c("wrmeanB"="Variant", "wrmeanW"="Non-variant")) + 
		scale_color_manual(values=c("wrmeanB"="red","wrmeanW"="blue","1"="red","0"="blue"), labels=c("wrmeanB"="Variant", "wrmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("wrmeanB"="red","wrmeanW"="blue"), labels=c("wrmeanB"="Variant", "wrmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=wr, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		labs(x="Mean clearance stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off


fig_wrmean_withprior <- shared_params_df %>% 
	# select(wrmeanB, wrmeanW) %>% 
	select(wrmeanB_trans, wrmeanW_trans) %>% 
	rename(wrmeanB=wrmeanB_trans, wrmeanW=wrmeanW_trans) %>%
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		# scale_color_manual(values=c("wrmeanB"="red","wrmeanW"="blue"), labels=c("wrmeanB"="Variant", "wrmeanW"="Non-variant")) + 
		scale_color_manual(values=c("wrmeanB"="red","wrmeanW"="blue","1"="red","0"="blue"), labels=c("wrmeanB"="Variant", "wrmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
		scale_fill_manual(values=c("wrmeanB"="red","wrmeanW"="blue"), labels=c("wrmeanB"="Variant", "wrmeanW"="Non-variant")) + 
		geom_jitter(data=meanvalsindiv, aes(x=wr, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$wrmean_prior, sd=prior_pars$wrsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		labs(x="Mean clearance stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

if(run_pars$trapfit==TRUE){
	fig_infdurmean <- shared_params_df %>% 
		mutate(infdurmeanB=wpmeanB_trans+wxmeanB_trans+wrmeanB_trans) %>%
		mutate(infdurmeanW=wpmeanW_trans+wxmeanW_trans+wrmeanW_trans) %>%
		select(infdurmeanB, infdurmeanW) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
			scale_x_continuous(limits=c(0,NA)) + 
			scale_color_manual(values=c("infdurmeanB"="red","infdurmeanW"="blue","1"="red","0"="blue"), labels=c("infdurmeanB"="Variant", "infdurmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
			scale_fill_manual(values=c("infdurmeanB"="red","infdurmeanW"="blue"), labels=c("infdurmeanB"="Variant", "infdurmeanW"="Non-variant")) + 
			geom_jitter(data=meanvalsindiv, aes(x=infdur, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
			labs(x="Mean acute infection duration (days)", y="Density") + 
			theme_minimal() + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off	
	} else {
	fig_infdurmean <- shared_params_df %>% 
		mutate(infdurmeanB=wpmeanB_trans+wrmeanB_trans) %>%
		mutate(infdurmeanW=wpmeanW_trans+wrmeanW_trans) %>%
		select(infdurmeanB, infdurmeanW) %>% 
		pivot_longer(everything()) %>% 
		ggplot(aes(x=value)) + 
			# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
			geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
			scale_x_continuous(limits=c(0,NA)) + 
			# scale_color_manual(values=c("infdurmeanB"="red","infdurmeanW"="blue"), labels=c("infdurmeanB"="Variant", "infdurmeanW"="Non-variant")) + 
			scale_color_manual(values=c("infdurmeanB"="red","infdurmeanW"="blue","1"="red","0"="blue"), labels=c("infdurmeanB"="Variant", "infdurmeanW"="Non-variant","1"="Variant","0"="Non-variant")) + 
			scale_fill_manual(values=c("infdurmeanB"="red","infdurmeanW"="blue"), labels=c("infdurmeanB"="Variant", "infdurmeanW"="Non-variant")) + 
			geom_jitter(data=meanvalsindiv, aes(x=infdur, y=jitterfactor*special, col=factor(special)), alpha=0.4, width=0, height=0.003)  + 
			# geom_point(data=meanvalsindiv, aes(x=wp+wr, y=jitterfactor*special, col=factor(special)))  + 
			labs(x="Mean acute infection duration (days)", y="Density") + 
			theme_minimal() + 
			# theme(legend.title=element_blank(), text=element_text(size=18)) + 
			theme(legend.position="none", text=element_text(size=18)) + 
			y_ticks_off + 
			grid_off	
	}

if(run_pars$trapfit==TRUE){
	fig_ct_trajectory_inference <- make_sample_trajectory_special_trap(shared_params_df, global_pars)
	} else {
	fig_ct_trajectory_inference <- make_sample_trajectory_special(shared_params_df, global_pars)		
	}

if(run_pars$trapfit==TRUE){
	fig_ct_trajectory_inference_green <- make_sample_trajectory_special_trap(shared_params_df, global_pars, specialfill="green")
	} else {
	fig_ct_trajectory_inference_green <- make_sample_trajectory_special(shared_params_df, global_pars, specialfill="green")		
	}

if(run_pars$trapfit==TRUE){
	fig_ct_trajectory_inference_black <- make_sample_trajectory_special_trap(shared_params_df, global_pars, specialfill="black")
	} else {
	fig_ct_trajectory_inference_black <- make_sample_trajectory_special(shared_params_df, global_pars, specialfill="black")		
	}

# fig_ct_trajectory_inference_geml <- make_sample_trajectory_special(shared_params_df, global_pars, ge=TRUE)



# fig_ct_trajectory_inference_rawspecial <- make_sample_trajectory_wtonly(shared_params_df, global_pars) + 
# 	geom_segment(data=filter(meanvalsindiv, special==1), aes(x=-wp, y=global_pars[["lod"]], xend=0, yend=global_pars[["lod"]]-dp), col="red", alpha=0.2) + 
# 	geom_segment(data=filter(meanvalsindiv, special==1), aes(x=0, y=global_pars[["lod"]]-dp, xend=wr, yend=global_pars[["lod"]]), col="red", alpha=0.2) +
# 	geom_point(data=filter(indiv_data, special==1), aes(x=t, y=y), col="red", alpha=0.2, size=1)

fig_ct_trajectory_inference + 
	# geom_segment(data=filter(meanvalsindiv, special==1), aes(x=-wp, y=global_pars[["lod"]], xend=0, yend=global_pars[["lod"]]-dp), col="red", alpha=0.2) + 
	# geom_segment(data=filter(meanvalsindiv, special==1), aes(x=0, y=global_pars[["lod"]]-dp, xend=wr, yend=global_pars[["lod"]]), col="red", alpha=0.2) +
	geom_point(data=filter(indiv_data, special==1), aes(x=t, y=y), col="red", alpha=0.2, size=1)


# fig_raw_trajectories <- with(as.list(global_pars),{
# meanvalsindiv %>% 
# 	arrange(special) %>%
# 	ggplot() + 
# 		geom_segment(aes(x=-Inf, y=lod, xend=-wp, yend=lod, col=factor(special), size=factor(special)), alpha=0.5) + 
# 		geom_segment(aes(x=-wp, y=lod, xend=0, yend=lod-dp, col=factor(special), size=factor(special)), alpha=0.5) + 
# 		geom_segment(aes(x=0, y=lod-dp, xend=wr, yend=lod, col=factor(special), size=factor(special)), alpha=0.5) + 
# 		geom_segment(aes(x=wr, y=lod, xend=Inf, yend=lod, col=factor(special), size=factor(special)), alpha=0.5) + 
# 		scale_color_manual(values=c("1"="red","0"="blue","Yes"="red","No"="blue"), labels=c("1"="Variant","0"="Non-variant","Yes"="Variant","No"="Non-variant")) + 
# 		scale_size_manual(values=c("1"=1, "0"=0.1), labels=c("1"="Variant","0"="Non-variant")) + 
# 		theme_minimal() + 
# 		theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), text=element_text(size=18)) + 
# 		labs(x="Time since min Ct (days)", y="Ct") + 
# 		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml"))
# })


fig_indiv_trajectories <- ggplot() + 
	geom_segment(
		data=filter(meanvalsindiv, special==0),
		aes(x=-wp, xend=0, y=40, yend=40-dp),
		size=0.4, alpha=0.2, col="blue") + 
	geom_segment(
		data=filter(meanvalsindiv, special==0),
		aes(x=0, xend=wr, y=40-dp, yend=40),
		size=0.4, alpha=0.2, col="blue") + 
	geom_segment(
		data=filter(meanvalsindiv, special==1),
		aes(x=-wp, xend=0, y=40, yend=40-dp),
		size=0.6, alpha=0.5, col="red") + 
	geom_segment(
		data=filter(meanvalsindiv, special==1),
		aes(x=0, xend=wr, y=40-dp, yend=40),
		size=0.6, alpha=0.5, col="red") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) + 
	theme_minimal() + 
	theme(text=element_text(size=14)) + 
	labs(x="Time (days)", y="Ct")

fig_indiv_trajectories_green <- ggplot() + 
	geom_segment(
		data=filter(meanvalsindiv, special==0),
		aes(x=-wp, xend=0, y=40, yend=40-dp),
		size=0.4, alpha=0.2, col="blue") + 
	geom_segment(
		data=filter(meanvalsindiv, special==0),
		aes(x=0, xend=wr, y=40-dp, yend=40),
		size=0.4, alpha=0.2, col="blue") + 
	geom_segment(
		data=filter(meanvalsindiv, special==1),
		aes(x=-wp, xend=0, y=40, yend=40-dp),
		size=0.6, alpha=0.5, col="green") + 
	geom_segment(
		data=filter(meanvalsindiv, special==1),
		aes(x=0, xend=wr, y=40-dp, yend=40),
		size=0.6, alpha=0.5, col="green") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) + 
	theme_minimal() + 
	theme(text=element_text(size=14)) + 
	labs(x="Time (days)", y="Ct")


