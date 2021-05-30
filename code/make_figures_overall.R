# For plotting output from two runs together (e.g. to compare a B.1.1.7 and B.1.429 fit); this one is a bit less refined than the rest of the code, since it relies on the ordering specified in set_run_pars, where the B.1.1.7 fits (vs non-VOI/VOC) are in slots 2, 6, and 10, and the B.1.429 fits (vs non-VOI/VOC) are in slots 4, 8, and 12. 

library(tidyverse) 
library(ggbeeswarm)

source('code/utils.R')
source('code/utils_private.R')

meanvalsindiv_combined <- tibble()
shared_params_df_combined <- tibble() 

for(run_pars_index in c(2,4)){ # c(6,8) c(10,12)

	# Import necessary information: -------------------------------------------
	run_pars <- run_pars_list[[run_pars_index]]

	savedir <- paste0("figures/run_pars_",run_pars_index,"/")
	savedir_parent <- paste0("figures/overall/")
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

	if(run_pars_index%in%c(2,6,10)){
		shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"B.1.1.7", TRUE~lineage)) 
	} else {
		shared_params_df_summary <- shared_params_df_summary %>% 
			mutate(lineage=case_when(lineage=="BX"~"B.1.429", TRUE~lineage)) 
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
			(run_pars_index%in%c(2,6,10) & special==1)~"B.1.1.7", 
			(run_pars_index%in%c(4,8,12) & special==1)~"B.1.429",
			TRUE~"NV"))
	meanvalsindiv_combined <- bind_rows(meanvalsindiv_combined, meanvalsindiv)

	meanvalsindiv_combined <- meanvalsindiv_combined %>% 
		group_by(id) %>% 
		summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr), infdur=mean(infdur), special=first(special), lineage=first(lineage))

	# Generate new trajectory figures with updated colors: --------------------
	if(run_pars$trapfit==TRUE){
		fig_ct_trajectory_inference <- make_sample_trajectory_special_trap(shared_params_df, global_pars, siglevel=0.95, specialfill=ifelse(run_pars_index%in%c(2,6,10),"red","green"))
		} else {
		fig_ct_trajectory_inference <- make_sample_trajectory_special(shared_params_df, global_pars, siglevel=0.95, specialfill=ifelse(run_pars_index%in%c(2,6,10),"red","green"))		
		}

	if(run_pars_index %in% c(2,6,10)){
		fig_ct_trajectory_b117 <- fig_ct_trajectory_inference
	} else {
		fig_ct_trajectory_b1429 <- fig_ct_trajectory_inference
	}

}

shared_params_df_summary <- shared_params_df_combined %>% 
	group_by(statistic, lineage) %>% 
	summarise(
		mean=mean(value),
		lwr=quantile(value,0.025),
		upr=quantile(value,0.975))
	
# -----------------------------------------------------------------------------
# Get ready for plotting (add x value columns first)
meanvalsindiv_combined_forplotting <- meanvalsindiv_combined %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="B.1.1.7"~2, lineage=="B.1.429"~3)) %>% 
	select(dp, wp, wr, infdur, lineage, xval) %>% 
	rename(id=infdur) %>% 
	pivot_longer(c("dp","wp","wr","id")) %>% 
	rename(statistic=name) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","B.1.1.7","B.1.429")))

shared_params_df_summary_forplotting <- shared_params_df_summary %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="B.1.1.7"~2, lineage=="B.1.429"~3)) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","B.1.1.7","B.1.429")))

# Set plotting parameters: 
whiskerwidth <- 0.075

# Splitting out: 
fig_whiskers_dp <- ggplot() + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
		aes(x=xval, y=global_pars[["lod"]]-value, col=lineage), alpha=0.4, size=2, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=global_pars[["lod"]]-mean,yend=global_pars[["lod"]]-mean), col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-lwr), col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-upr,yend=global_pars[["lod"]]-upr), col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval,xend=xval,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-upr), col="gray20") + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/VOC","B.1.1.7","B.1.429")) + 
	scale_color_manual(values=c("blue","red","green")) + 
	labs(title="Peak viral concentration", x=element_blank(), y="Ct") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14), legend.position="none") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) 

fig_whiskers_wp <- ggplot() + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
		aes(x=xval, y=value, col=lineage), alpha=0.4, size=2, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20") + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/VOC","B.1.1.7","B.1.429")) + 
	scale_color_manual(values=c("blue","red","green")) + 
	labs(title="Proliferation time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14), legend.position="none")

fig_whiskers_wr <- ggplot() + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
		aes(x=xval, y=value, col=lineage), alpha=0.4, size=2, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20") + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/VOC","B.1.1.7","B.1.429")) + 
	scale_color_manual(values=c("blue","red","green")) + 
	labs(title="Clearance time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14), legend.position="none")

fig_whiskers_infdur <- ggplot() + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
		aes(x=xval, y=value, col=lineage), alpha=0.4, size=2, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20") +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20") + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/VOC","B.1.1.7","B.1.429")) + 
	scale_color_manual(values=c("blue","red","green")) + 
	labs(title="Acute infection duration", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14), legend.position="none")


# Save: -----------------------------------------------------------------------
ggsave(fig_whiskers_dp, file=paste0(savedir_parent,"whiskers_dp.pdf"), width=4, height=4)
ggsave(fig_whiskers_wp, file=paste0(savedir_parent,"whiskers_wp.pdf"), width=4, height=4)
ggsave(fig_whiskers_wr, file=paste0(savedir_parent,"whiskers_wr.pdf"), width=4, height=4)
ggsave(fig_whiskers_infdur, file=paste0(savedir_parent,"whiskers_infdur.pdf"), width=4, height=4)
ggsave(fig_ct_trajectory_b117, file=paste0(savedir_parent,"ct_trajectory_b117.pdf"), width=8, height=5)
ggsave(fig_ct_trajectory_b1429, file=paste0(savedir_parent,"ct_trajectory_b1429.pdf"), width=8, height=5)

# Some KS tests: --------------------------------------------------------------
# B117:
ks_b117_dp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.1.7")$dp,
	filter(meanvalsindiv_combined,lineage=="NV")$dp
	)$p.value

ks_b117_wp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.1.7")$wp,
	filter(meanvalsindiv_combined,lineage=="NV")$wp
	)$p.value

ks_b117_wr <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.1.7")$wr,
	filter(meanvalsindiv_combined,lineage=="NV")$wr
	)$p.value

ks_b117_infdur <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.1.7")$infdur,
	filter(meanvalsindiv_combined,lineage=="NV")$infdur
	)$p.value

c(ks_b117_dp, ks_b117_wp, ks_b117_wr, ks_b117_infdur)

# B1429:
ks_b1429_dp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.429")$dp,
	filter(meanvalsindiv_combined,lineage=="NV")$dp
	)$p.value

ks_b1429_wp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.429")$wp,
	filter(meanvalsindiv_combined,lineage=="NV")$wp
	)$p.value

ks_b1429_wr <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.429")$wr,
	filter(meanvalsindiv_combined,lineage=="NV")$wr
	)$p.value

ks_b1429_infdur <- ks.test(
	filter(meanvalsindiv_combined,lineage=="B.1.429")$infdur,
	filter(meanvalsindiv_combined,lineage=="NV")$infdur
	)$p.value

c(ks_b1429_dp, ks_b1429_wp, ks_b1429_wr, ks_b1429_infdur)

# Summarize individual values: ------------------------------------------------

indiv_table_ct <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		ct_min=min(global_pars[["lod"]]-dp),
		ct_lwr=quantile(global_pars[["lod"]]-dp, 0.25),
		ct_median=median(global_pars[["lod"]]-dp),
		ct_upr=quantile(global_pars[["lod"]]-dp, 0.75),
		ct_max=max(global_pars[["lod"]]-dp)
		)

indiv_table_geml <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		geml_min=min(convert_Ct_logGEML(global_pars[["lod"]]-dp)),
		geml_lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dp), 0.25),
		geml_median=median(convert_Ct_logGEML(global_pars[["lod"]]-dp)),
		geml_upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dp), 0.75),
		geml_max=max(convert_Ct_logGEML(global_pars[["lod"]]-dp))
		)

indiv_table_wp <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		wp_min=min(wp),
		wp_lwr=quantile(wp, 0.25),
		wp_median=median(wp),
		wp_upr=quantile(wp, 0.75),
		wp_max=max(wp)
		)

indiv_table_wr <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		wr_min=min(wr),
		wr_lwr=quantile(wr, 0.25),
		wr_median=median(wr),
		wr_upr=quantile(wr, 0.75),
		wr_max=max(wr)
		)

indiv_table_infdur <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		infdur_min=min(infdur),
		infdur_lwr=quantile(infdur, 0.25),
		infdur_median=median(infdur),
		infdur_upr=quantile(infdur, 0.75),
		infdur_max=max(infdur)
		)

print(indiv_table_ct)
print(indiv_table_geml)
print(indiv_table_wp)
print(indiv_table_wr)
print(indiv_table_infdur)
