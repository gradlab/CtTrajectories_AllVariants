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


# runsets <- c(23) # uninformaative
runsets <- c(30) # informed
# runsets <- c(31) # low
for(run_pars_index in runsets){ 

	# Import necessary information: -------------------------------------------
	run_pars <- run_pars_list[[run_pars_index]]

	savedir <- paste0("figures/run_pars_",run_pars_index,"/")

	if(run_pars_index==23){
		savedir_parent <- paste0("figures/overall/uninformative/withvax/")
	}
	else if(run_pars_index==30){
		savedir_parent <- paste0("figures/overall/informed/withvax/")
	} else if(run_pars_index==31){
		savedir_parent <- paste0("figures/overall/low/withvax/")
	}
	
	savedir_output <- paste0("output/run_pars_",run_pars_index,"/")
	
	source('code/fit_posteriors_preamble.R')
	load(file=paste0(extdrive,savedir_output,"ct_fit.RData"))
	source('code/extract_parameters.R')


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

	shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"Vaccinated", TRUE~"Unvaccinated")) 

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
			(special==1)~"Vaccinated",
			(special==0)~"Unvaccinated"))
	meanvalsindiv_combined <- bind_rows(meanvalsindiv_combined, meanvalsindiv)

	meanvalsindiv_combined <- meanvalsindiv_combined %>% 
		group_by(id,lineage) %>% 
		summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr), infdur=mean(infdur), special=first(special)) %>%
		left_join(breakthroughdf, by="id")

	params_df <- params_df %>% 
		mutate(lineage=case_when(
			(special==1)~"Vaccinated",
			(special==0)~"Unvaccinated"))
	params_df_combined <- bind_rows(params_df_combined, params_df)

}

fig_ct_trajectory_deltavax <- make_sample_trajectory_special(shared_params_df, global_pars, siglevel=0.95, specialfill="red", referencefill="blue")	

shared_params_df_summary <- shared_params_df_combined %>% 
	group_by(statistic, lineage) %>% 
	summarise(
		mean=mean(value),
		lwr=quantile(value,0.025),
		upr=quantile(value,0.975))


# -----------------------------------------------------------------------------
# Make some plots with just the variants (no epsilon): 
meanvalsindiv_combined_forplotting_deltavax <- meanvalsindiv_combined %>% 
	ungroup() %>% 
	mutate(xval=case_when(lineage=="Unvaccinated"~1, lineage=="Vaccinated"~2)) %>% 
	select(dp, wp, wr, infdur, lineage, xval, breakthrough) %>% 
	rename(id=infdur) %>% 
	pivot_longer(c("dp","wp","wr","id")) %>% 
	rename(statistic=name) %>% 
	mutate(lineage=factor(lineage, levels=c("Unvaccinated","Vaccinated")))

shared_params_df_summary_forplotting_deltavax <- shared_params_df_summary %>% 
	mutate(xval=case_when(lineage=="Unvaccinated"~1, lineage=="Vaccinated"~2)) %>% 
	mutate(lineage=factor(lineage, levels=c("Unvaccinated","Vaccinated")))

# Set plotting parameters: 
whiskerwidth <- 0.05
jitterwidth <- 0.05

fig_whiskers_dp_deltavax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
	# 	aes(x=xval, y=global_pars[["lod"]]-value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_deltavax, statistic=="dp"),
		aes(x=xval, y=global_pars[["lod"]]-value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=4) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="dp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=global_pars[["lod"]]-mean,yend=global_pars[["lod"]]-mean), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-lwr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-upr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="dp"), 
		aes(x=xval,xend=xval,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("blue","red")) + 
	labs(title="Peak viral concentration", x=element_blank(), y="Ct") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) 

fig_whiskers_wp_deltavax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_deltavax, statistic=="wp"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=4) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wp"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("blue","red")) + 
	labs(title="Proliferation time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_wr_deltavax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_deltavax, statistic=="wr"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=4) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wr"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="wr"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("blue","red")) + 
	labs(title="Clearance time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_infdur_deltavax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_deltavax, statistic=="id"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=4) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="id"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_deltavax, statistic=="id"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("blue","red")) + 
	labs(title="Acute infection duration", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

ggsave(fig_whiskers_dp_deltavax, file=paste0(savedir_parent,"whiskers_dp_deltavax.pdf"), width=6, height=4)
ggsave(fig_whiskers_wp_deltavax, file=paste0(savedir_parent,"whiskers_wp_deltavax.pdf"), width=6, height=4)
ggsave(fig_whiskers_wr_deltavax, file=paste0(savedir_parent,"whiskers_wr_deltavax.pdf"), width=6, height=4)
ggsave(fig_whiskers_infdur_deltavax, file=paste0(savedir_parent,"whiskers_infdur_deltavax.pdf"), width=6, height=4)
ggsave(fig_ct_trajectory_deltavax, file=paste0(savedir_parent,"ct_trajectory_deltavax.pdf"), width=8, height=5)


ggsave(fig_whiskers_dp_deltavax, file=paste0(savedir_parent,"whiskers_dp_deltavax.png"), width=6, height=4)
ggsave(fig_whiskers_wp_deltavax, file=paste0(savedir_parent,"whiskers_wp_deltavax.png"), width=6, height=4)
ggsave(fig_whiskers_wr_deltavax, file=paste0(savedir_parent,"whiskers_wr_deltavax.png"), width=6, height=4)
ggsave(fig_whiskers_infdur_deltavax, file=paste0(savedir_parent,"whiskers_infdur_deltavax.png"), width=6, height=4)
ggsave(fig_ct_trajectory_deltavax, file=paste0(savedir_parent,"ct_trajectory_deltavax.png"), width=8, height=5)