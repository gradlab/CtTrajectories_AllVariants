library(tidyverse) 
library(ggbeeswarm)

source('code/utils.R')
source('code/utils_private.R')

breakthroughdf <- ct_dat_refined %>% 
	select(id=PersonID, breakthrough=VaccineBreakthrough) %>%
	group_by(id) %>%
	slice(1)

symptomdf <- ct_dat_refined %>% 
	select(id=PersonID, symptoms=Symptoms) %>%
	group_by(id) %>%
	slice(1)


meanvalsindiv_combined <- tibble()
params_df_combined <- tibble()
shared_params_df_combined <- tibble() 

# c(2,4,6): uninformative
# c(8,10,12): informed
# c(14,16,18): low

# runsets <- c(2,4,6,24)
runsets <- c(8,10,12,28)
# runsets <- c(14,16,18,29)
for(run_pars_index in runsets){ 

	# Import necessary information: -------------------------------------------
	run_pars <- run_pars_list[[run_pars_index]]

	savedir <- paste0("figures/run_pars_",run_pars_index,"/")

	if(run_pars_index %in% c(2,4,6)){
		savedir_parent <- paste0("figures/overall/uninformative/withvax/")
	} else if(run_pars_index %in% c(8,10,12)){
		savedir_parent <- paste0("figures/overall/informed/withvax/")	
	} else if(run_pars_index %in% c(14,16,18)){
		savedir_parent <- paste0("figures/overall/low/withvax/")	
	}
	
	savedir_output <- paste0("output/run_pars_",run_pars_index,"/")
	
	source('code/fit_posteriors_preamble.R')
	load(file=paste0(extdrive,savedir_output,"ct_fit.RData"))
	source('code/extract_parameters.R')

	# Plot the Ct trajectories individually: ----------------------------------
	if(run_pars_index %in% c(2,8,14)){

		# fig_ct_fit_nonVOIVOC <- plot_ct_fit(filter(params_df,special==0), global_pars, filter(indiv_data,special==0), ctalpha=0.01, specialcol="cornflowerblue", xlim=c(-20,30), nvlabel="Non-VOI/VOC", ntraces=100)	+ 
		# 	geom_rect(data=filter(breakthroughdf, id%in%(filter(params_df,special==0)$id)),aes(fill=breakthrough),xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, col="white", alpha=0.2) + 
		# 	scale_fill_manual(values=c("Yes"="Gray","No"="White"), guide="none")

		figlist_ct_fit_nonVOIVOC <- plot_ct_fit_split(filter(params_df,special==0), global_pars, filter(indiv_data,special==0), ctalpha=0.01, specialcol="cornflowerblue", xlim=c(-20,30), nvlabel="Non-VOI/VOC", nsplits=2, shade_ids=filter(breakthroughdf, breakthrough=="Yes")$id, ntraces=100)

		fig_ct_fit_alpha <- plot_ct_fit(filter(params_df,special==1), global_pars, filter(indiv_data,special==1), ctalpha=0.01, specialcol="red2", xlim=c(-20,30), vlabel="Alpha", ntraces=100)	+ 
			geom_rect(data=filter(breakthroughdf, id%in%(filter(params_df,special==1)$id)),aes(fill=breakthrough),xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, col="white", alpha=0.2) + 
			scale_fill_manual(values=c("Yes"="Gray","No"="White"), guide="none")

	} else if(run_pars_index %in% c(4,10,16)){

		fig_ct_fit_epsilon <- plot_ct_fit(filter(params_df,special==1), global_pars, filter(indiv_data,special==1), ctalpha=0.01, specialcol="forestgreen", xlim=c(-20,30), vlabel="Epsilon", ntraces=100)	+ 
			geom_rect(data=filter(breakthroughdf, id%in%(filter(params_df,special==1)$id)),aes(fill=breakthrough),xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, col="white", alpha=0.2) + 
			scale_fill_manual(values=c("Yes"="Gray","No"="White"), guide="none")

	} else if(run_pars_index %in% c(6,12,18)){

		fig_ct_fit_delta <- plot_ct_fit(filter(params_df,special==1), global_pars, filter(indiv_data,special==1), ctalpha=0.01, specialcol="orange", xlim=c(-20,30), vlabel="Delta", ntraces=100)	+ 
			geom_rect(data=filter(breakthroughdf, id%in%(filter(params_df,special==1)$id)),aes(fill=breakthrough),xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, col="white", alpha=0.2) + 
			scale_fill_manual(values=c("Yes"="Gray","No"="White"), guide="none")

	} else if(run_pars_index %in% c(24,28,29)){

		figlist_ct_fit_unvaccinated <- plot_ct_fit_split(filter(params_df,special==0), global_pars, filter(indiv_data,special==0), ctalpha=0.01, specialcol="darkblue", basecol="black", xlim=c(-20,30), nvlabel="Unvaccinated", nsplits=4, shade_ids=filter(breakthroughdf, breakthrough=="Yes")$id, ntraces=100)

		figlist_ct_fit_vaccinated <- plot_ct_fit_split(filter(params_df,special==1), global_pars, filter(indiv_data,special==1), ctalpha=0.01, specialcol="darkblue", basecol="black", xlim=c(-20,30), vlabel="Vaccinated", nsplits=1, shade_ids=filter(breakthroughdf, breakthrough=="Yes")$id, ntraces=100)

	}
	

	# Amend params_df: --------------------------------------------------------

	shared_params_df_summary <- shared_params_df %>% 
		select(dpmeanW_trans, dpmeanB_trans, wpmeanW_trans, wpmeanB_trans, wrmeanW_trans, wrmeanB_trans) %>% 
		mutate(idmeanW_trans=wpmeanW_trans+wrmeanW_trans) %>% 
		mutate(idmeanB_trans=wpmeanB_trans+wrmeanB_trans) %>% 
		pivot_longer(everything()) %>% 
		mutate(statistic=substr(name,1,2)) %>% 
		mutate(lineage=substr(name,7,7)) %>% 
		mutate(lineage=case_when(lineage=="B"~"BX", TRUE~"NV")) %>% 
		select(statistic, lineage, value)

	if(run_pars_index%in%c(2,8,14)){
		shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"Alpha", TRUE~lineage)) 
	} else if(run_pars_index%in%c(4,10,16)){
		shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"Epsilon", TRUE~lineage)) 
	} else if(run_pars_index%in%c(6,12,18)){
		shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"Delta", TRUE~lineage)) 
	} else if(run_pars_index%in%c(24,28,29)){
		shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"Vaccinated", TRUE~"Unvaccinated")) 
	}

	shared_params_df_combined <- bind_rows(shared_params_df_combined, shared_params_df_summary)

	# Store mean individual values: -------------------------------------------
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
	meanvalsindiv <- meanvalsindiv %>% 
		mutate(lineage=case_when(
			(run_pars_index%in%c(2,8,14) & special==1)~"Alpha", 
			(run_pars_index%in%c(4,10,16) & special==1)~"Epsilon",
			(run_pars_index%in%c(6,12,18) & special==1)~"Delta",
			(run_pars_index%in%c(24,28,29) & special==1)~"Vaccinated",
			(run_pars_index%in%c(24,28,29) & special==0)~"Unvaccinated",
			TRUE~"NV"))
	meanvalsindiv_combined <- bind_rows(meanvalsindiv_combined, meanvalsindiv)

	meanvalsindiv_combined <- meanvalsindiv_combined %>% 
		group_by(id,lineage) %>% 
		summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr), infdur=mean(infdur), special=first(special)) %>%
		left_join(breakthroughdf, by="id")

	params_df <- params_df %>% 
		mutate(lineage=case_when(
			(run_pars_index%in%c(2,8,14) & special==1)~"Alpha", 
			(run_pars_index%in%c(4,10,16) & special==1)~"Epsilon",
			(run_pars_index%in%c(6,12,18) & special==1)~"Delta",
			(run_pars_index%in%c(24,28,29) & special==1)~"Vaccinated",
			(run_pars_index%in%c(24,28,29) & special==0)~"Unvaccinated",
			TRUE~"NV"))
	params_df_combined <- bind_rows(params_df_combined, params_df)

	# Generate new trajectory figures with updated colors: --------------------
	if(run_pars_index %in% c(2,8,14)){
		specialfillval <- "red2"
		referencefillval <- "cornflowerblue"
	} else if(run_pars_index %in% c(4,10,16)){
		specialfillval <- "forestgreen"
		referencefillval <- "cornflowerblue"
	} else if(run_pars_index %in% c(6,12,18)){
		specialfillval <- "orange"
		referencefillval <- "cornflowerblue"
	} else if(run_pars_index %in% c(24,28,29)){
		specialfillval <- "darkblue"
		referencefillval <- "forestgreen"
	}

# blue -> cornflowerblue
# red -> red2
# purple -> orange
# green -> forestgreen
# maroon -> darkblue

	if(run_pars$trapfit==TRUE){
		fig_ct_trajectory_inference <- make_sample_trajectory_special_trap(shared_params_df, global_pars, siglevel=0.95, specialfill=specialfillval, referencefill=referencefillval)
		} else {
		fig_ct_trajectory_inference <- make_sample_trajectory_special(shared_params_df, global_pars, siglevel=0.95, specialfill=specialfillval, referencefill=referencefillval)		
		}

	if(run_pars_index %in% c(2,8,14)){
		fig_ct_trajectory_alpha <- fig_ct_trajectory_inference
	} else if(run_pars_index %in% c(4,10,16)){
		fig_ct_trajectory_epsilon <- fig_ct_trajectory_inference
	} else if(run_pars_index %in% c(6,12,18)){
		fig_ct_trajectory_delta <- fig_ct_trajectory_inference
	} else if(run_pars_index %in% c(24,28,29)){
		fig_ct_trajectory_vax <- fig_ct_trajectory_inference
	}

}

shared_params_df_summary <- shared_params_df_combined %>% 
	group_by(statistic, lineage) %>% 
	summarise(
		mean=mean(value),
		lwr=quantile(value,0.025),
		upr=quantile(value,0.975))

# system("say Done with making figures overall")