
# -----------------------------------------------------------------------------
# Get ready for plotting (add x value columns first)
meanvalsindiv_combined_forplotting <- meanvalsindiv_combined %>% 
	ungroup() %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="Alpha"~2, lineage=="Epsilon"~3, lineage=="Delta"~4, lineage=="Unvaccinated"~5, lineage=="Vaccinated"~6)) %>% 
	select(dp, wp, wr, infdur, lineage, xval, breakthrough) %>% 
	rename(id=infdur) %>% 
	pivot_longer(c("dp","wp","wr","id")) %>% 
	rename(statistic=name) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated")))

shared_params_df_summary_forplotting <- shared_params_df_summary %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="Alpha"~2, lineage=="Epsilon"~3, lineage=="Delta"~4, lineage=="Unvaccinated"~5, lineage=="Vaccinated"~6)) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated")))

# Set plotting parameters: 
whiskerwidth <- 0.075
jitterwidth <- 0.05

# an attempt at all together: 
fig_whiskers_all <- ggplot() + 
	geom_jitter(data=meanvalsindiv_combined_forplotting,
		aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4) + 
	geom_segment(data=shared_params_df_summary_forplotting, 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean,col=lineage)) +
	geom_segment(data=shared_params_df_summary_forplotting, 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=lwr,yend=lwr,col=lineage)) +
	geom_segment(data=shared_params_df_summary_forplotting, 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=upr,yend=upr,col=lineage)) +
	geom_segment(data=shared_params_df_summary_forplotting, 
		aes(x=xval,xend=xval,y=lwr,yend=upr,col=lineage)) +
	facet_wrap(~statistic, nrow=1, scales="free",
		strip.position = "left", labeller = as_labeller(c(dp="DP", wp="WP", wr="WR", id="INFDUR"))) + 
	theme_minimal() + 
	theme(strip.background=element_blank(), strip.placement="outside")


meanvalsindiv_combined %>% 
	ggplot(aes(x=global_pars[["lod"]]-dp)) + 
		geom_histogram(aes(y=..density..), binwidth=1, col="black", fill="lightgray") + 
		theme_minimal() + 
		facet_wrap(~lineage)

# Splitting out: 
fig_whiskers_dp <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
	# 	aes(x=xval, y=global_pars[["lod"]]-value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
		aes(x=xval, y=global_pars[["lod"]]-value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=global_pars[["lod"]]-mean,yend=global_pars[["lod"]]-mean), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-lwr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-upr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="dp"), 
		aes(x=xval,xend=xval,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:6, minor_breaks=1:6, label=c("Non-VOI/\nVOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","forestgreen","orange","black","darkblue")) + 
	labs(title="Peak viral concentration", x=element_blank(), y="Ct") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) 

fig_whiskers_wp <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wp"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:6, minor_breaks=1:6, label=c("Non-VOI/\nVOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","forestgreen","orange","black","darkblue")) + 
	labs(title="Proliferation time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_wr <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="wr"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:6, minor_breaks=1:6, label=c("Non-VOI/\nVOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","forestgreen","orange","black","darkblue")) + 
	labs(title="Clearance time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_infdur <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting, statistic=="id"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:6, minor_breaks=1:6, label=c("Non-VOI/\nVOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","forestgreen","orange","black","darkblue")) + 
	labs(title="Acute infection duration", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")


# -----------------------------------------------------------------------------
# Make some plots without epsilon: 
meanvalsindiv_combined_forplotting_noepsilon <- meanvalsindiv_combined %>% 
	ungroup() %>% 
	filter(lineage!="Epsilon") %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="Alpha"~2, lineage=="Delta"~3, lineage=="Unvaccinated"~4, lineage=="Vaccinated"~5)) %>% 
	select(dp, wp, wr, infdur, lineage, xval, breakthrough) %>% 
	rename(id=infdur) %>% 
	pivot_longer(c("dp","wp","wr","id")) %>% 
	rename(statistic=name) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","Alpha","Delta","Unvaccinated","Vaccinated")))

shared_params_df_summary_forplotting_noepsilon <- shared_params_df_summary %>% 
	filter(lineage!="Epsilon") %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="Alpha"~2, lineage=="Delta"~3, lineage=="Unvaccinated"~4, lineage=="Vaccinated"~5)) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","Alpha","Delta","Unvaccinated","Vaccinated")))

# Set plotting parameters: 
whiskerwidth <- 0.075
jitterwidth <- 0.05

fig_whiskers_dp_noepsilon <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
	# 	aes(x=xval, y=global_pars[["lod"]]-value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_noepsilon, statistic=="dp"),
		aes(x=xval, y=global_pars[["lod"]]-value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="dp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=global_pars[["lod"]]-mean,yend=global_pars[["lod"]]-mean), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-lwr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-upr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="dp"), 
		aes(x=xval,xend=xval,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:5, minor_breaks=1:5, label=c("Non-VOI/\nVOC","Alpha","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange","forestgreen","darkblue")) + 
	labs(title="Peak viral concentration", x=element_blank(), y="Ct") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) 

fig_whiskers_wp_noepsilon <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_noepsilon, statistic=="wp"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wp"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:5, minor_breaks=1:5, label=c("Non-VOI/\nVOC","Alpha","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange","forestgreen","darkblue")) + 
	labs(title="Proliferation time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_wr_noepsilon <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_noepsilon, statistic=="wr"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wr"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="wr"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:5, minor_breaks=1:5, label=c("Non-VOI/\nVOC","Alpha","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange","forestgreen","darkblue")) + 
	labs(title="Clearance time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_infdur_noepsilon <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=2) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_noepsilon, statistic=="id"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=2, cex=1.7) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="id"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_noepsilon, statistic=="id"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:5, minor_breaks=1:5, label=c("Non-VOI/\nVOC","Alpha","Delta","Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange","forestgreen","darkblue")) + 
	labs(title="Acute infection duration", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

# -----------------------------------------------------------------------------
# Make some plots with just the variants (no epsilon): 
meanvalsindiv_combined_forplotting_justvariants <- meanvalsindiv_combined %>% 
	ungroup() %>% 
	filter(lineage %in% c("NV","Alpha","Delta")) %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="Alpha"~2, lineage=="Delta"~3)) %>% 
	select(dp, wp, wr, infdur, lineage, xval, breakthrough) %>% 
	rename(id=infdur) %>% 
	pivot_longer(c("dp","wp","wr","id")) %>% 
	rename(statistic=name) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","Alpha","Delta")))

shared_params_df_summary_forplotting_justvariants <- shared_params_df_summary %>% 
	filter(lineage %in% c("NV","Alpha","Delta")) %>% 
	mutate(xval=case_when(lineage=="NV"~1, lineage=="Alpha"~2, lineage=="Delta"~3)) %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, levels=c("Non-VOI/VOC","Alpha","Delta")))

# Set plotting parameters: 
whiskerwidth <- 0.11
jitterwidth <- 0.05

fig_whiskers_dp_justvariants <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
	# 	aes(x=xval, y=global_pars[["lod"]]-value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvariants, statistic=="dp"),
		aes(x=xval, y=global_pars[["lod"]]-value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="dp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=global_pars[["lod"]]-mean,yend=global_pars[["lod"]]-mean), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-lwr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-upr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="dp"), 
		aes(x=xval,xend=xval,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/\nVOC","Alpha","Delta")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange")) + 
	labs(title="Peak viral concentration", x=element_blank(), y="Ct") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) 

fig_whiskers_wp_justvariants <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvariants, statistic=="wp"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wp"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/\nVOC","Alpha","Delta")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange")) + 
	labs(title="Proliferation time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_wr_justvariants <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvariants, statistic=="wr"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wr"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="wr"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/\nVOC","Alpha","Delta")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange")) + 
	labs(title="Clearance time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_infdur_justvariants <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvariants, statistic=="id"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="id"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvariants, statistic=="id"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:3, minor_breaks=1:3, label=c("Non-VOI/\nVOC","Alpha","Delta")) + 
	scale_color_manual(values=c("cornflowerblue","red2","orange")) + 
	labs(title="Acute infection duration", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

# -----------------------------------------------------------------------------
# Make some plots with just vaccine status (no variants: 
meanvalsindiv_combined_forplotting_justvax <- meanvalsindiv_combined %>% 
	ungroup() %>% 
	filter(lineage %in% c("Unvaccinated","Vaccinated")) %>% 
	mutate(xval=case_when(lineage=="Unvaccinated"~1, lineage=="Vaccinated"~2)) %>% 
	select(dp, wp, wr, infdur, lineage, xval, breakthrough) %>% 
	rename(id=infdur) %>% 
	pivot_longer(c("dp","wp","wr","id")) %>% 
	rename(statistic=name) %>% 
	mutate(lineage=factor(lineage, levels=c("Unvaccinated","Vaccinated")))

shared_params_df_summary_forplotting_justvax <- shared_params_df_summary %>% 
	filter(lineage %in% c("Unvaccinated","Vaccinated")) %>%
	mutate(xval=case_when(lineage=="Unvaccinated"~1, lineage=="Vaccinated"~2)) %>% 
	mutate(lineage=factor(lineage, levels=c("Unvaccinated","Vaccinated")))

# Set plotting parameters: 
whiskerwidth <- 0.075
jitterwidth <- 0.05

fig_whiskers_dp_justvax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="dp"),
	# 	aes(x=xval, y=global_pars[["lod"]]-value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvax, statistic=="dp"),
		aes(x=xval, y=global_pars[["lod"]]-value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3.5) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="dp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=global_pars[["lod"]]-mean,yend=global_pars[["lod"]]-mean), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-lwr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="dp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=global_pars[["lod"]]-upr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="dp"), 
		aes(x=xval,xend=xval,y=global_pars[["lod"]]-lwr,yend=global_pars[["lod"]]-upr), col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("forestgreen","darkblue")) + 
	labs(title="Peak viral concentration", x=element_blank(), y="Ct") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none") + 
	scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) 

fig_whiskers_wp_justvax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wp"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvax, statistic=="wp"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3.5) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wp"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wp"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wp"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("forestgreen","darkblue")) + 
	labs(title="Proliferation time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_wr_justvax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="wr"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvax, statistic=="wr"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3.5) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wr"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wr"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="wr"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("forestgreen","darkblue")) + 
	labs(title="Clearance time", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")

fig_whiskers_infdur_justvax <- ggplot() + 
	# geom_jitter(data=filter(meanvalsindiv_combined_forplotting, statistic=="id"),
	# 	aes(x=xval, y=value, col=lineage), width=jitterwidth, height=0, alpha=0.4, size=3) + 
	geom_beeswarm(data=filter(meanvalsindiv_combined_forplotting_justvax, statistic=="id"),
		aes(x=xval, y=value, col=lineage, shape=breakthrough), alpha=0.4, size=3, cex=3.5) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="id"), 
		aes(x=xval-whiskerwidth,xend=xval+whiskerwidth,y=mean,yend=mean),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=lwr,yend=lwr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="id"), 
		aes(x=xval-whiskerwidth/2,xend=xval+whiskerwidth/2,y=upr,yend=upr),col="gray20", size=.75) +
	geom_segment(data=filter(shared_params_df_summary_forplotting_justvax, statistic=="id"), 
		aes(x=xval,xend=xval,y=lwr,yend=upr),col="gray20", size=.75) + 
	scale_x_continuous(breaks=1:2, minor_breaks=1:2, label=c("Unvaccinated","Vaccinated")) + 
	scale_color_manual(values=c("forestgreen","darkblue")) + 
	labs(title="Acute infection duration", x=element_blank(), y="Days") + 
	theme_minimal() + 
	theme(text=element_text(size=14), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1.1), legend.position="none")



# Plot individual-level probabilities of Ct under some value: ------------------

punder_labs <- params_df_combined %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC", TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, c("Non-VOI/VOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>%
	mutate(MinCt=global_pars[["lod"]]-dp) %>% 
	mutate(MinCtUnder15=case_when(MinCt<15 ~ 1,TRUE~0)) %>% 
	group_by(lineage) %>% 
	summarise(punder=round(sum(MinCtUnder15)/n(),3)) %>% 
	mutate(punder=paste0("Proportion under\n15 Ct: ",as.character(punder)))

fig_lineage_hists <- params_df_combined %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC", TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, c("Non-VOI/VOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>% 
	ggplot(aes(x=global_pars[["lod"]]-dp)) + 
		geom_histogram(aes(y=..density..),binwidth=1, col="black", fill="lightgray") + 
		geom_density() + 
		geom_text(data=punder_labs, aes(label=punder), x=7.5, y=0.07) + 
		# theme_void() + 
		geom_vline(xintercept=15) + 
		facet_wrap(~lineage) + 
		theme(text=element_text(size=14)) + 
		labs(x="Ct",y="Density")

punder_labs_justvariants <- params_df_combined %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC", TRUE~lineage)) %>% 
	filter(lineage %in% c("Non-VOI/VOC", "Alpha","Delta")) %>% 
	mutate(lineage=factor(lineage, c("Non-VOI/VOC","Alpha","Delta"))) %>%
	mutate(MinCt=global_pars[["lod"]]-dp) %>% 
	mutate(MinCtUnder15=case_when(MinCt<15 ~ 1,TRUE~0)) %>% 
	group_by(lineage) %>% 
	summarise(punder=round(sum(MinCtUnder15)/n(),3)) %>% 
	mutate(punder=paste0("Proportion under\n15 Ct: ",as.character(punder)))

fig_lineage_hists_justvariants <- params_df_combined %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC", TRUE~lineage)) %>% 
	filter(lineage %in% c("Non-VOI/VOC", "Alpha","Delta")) %>% 
	mutate(lineage=factor(lineage, c("Non-VOI/VOC","Alpha","Delta"))) %>% 
	ggplot(aes(x=global_pars[["lod"]]-dp)) + 
		geom_histogram(aes(y=..density..),binwidth=1, col="gray", fill="white") + 
		geom_density(adjust=3) + 
		# geom_text(data=punder_labs_justvariants, aes(label=punder), x=7.5, y=0.07) + 
		theme_minimal() + 
		theme(panel.grid=element_blank()) + 
		geom_vline(xintercept=15, size=0.2, lty="dashed") + 
		facet_wrap(~lineage) + 
		theme(text=element_text(size=18)) + 
		labs(x="Peak Ct",y="Density") + 
		scale_x_reverse(limits=c(40,0)) 
		# scale_x_continuous(breaks=NULL)
		# scale_x_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(Peak~log[10]~RNA~copies/ml))) 

punder_labs_justvax <- params_df_combined %>% 
	filter(lineage %in% c("Vaccinated", "Unvaccinated")) %>% 
	mutate(lineage=factor(lineage, c("Vaccinated", "Unvaccinated"))) %>%
	mutate(MinCt=global_pars[["lod"]]-dp) %>% 
	mutate(MinCtUnder15=case_when(MinCt<15 ~ 1,TRUE~0)) %>% 
	group_by(lineage) %>% 
	summarise(punder=round(sum(MinCtUnder15)/n(),3)) %>% 
	mutate(punder=paste0("Proportion under\n15 Ct: ",as.character(punder)))

fig_lineage_hists_justvax <- params_df_combined %>% 
	filter(lineage %in% c("Vaccinated", "Unvaccinated")) %>% 
	mutate(lineage=factor(lineage, c("Vaccinated", "Unvaccinated"))) %>% 
	ggplot(aes(x=global_pars[["lod"]]-dp)) + 
		geom_histogram(aes(y=..density..),binwidth=1, col="gray", fill="white") + 
		geom_density() + 
		# geom_text(data=punder_labs_justvariants, aes(label=punder), x=7.5, y=0.07) + 
		theme_minimal() + 
		theme(panel.grid=element_blank()) + 
		geom_vline(xintercept=15, size=0.2, lty="dashed") + 
		facet_wrap(~lineage) + 
		theme(text=element_text(size=18)) + 
		labs(x="Peak Ct",y="Density") + 
		scale_x_reverse(limits=c(40,0)) 
		# scale_x_continuous(breaks=NULL)
		# scale_x_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(Peak~log[10]~RNA~copies/ml))) 

fig_lineage_densities <- params_df_combined %>% 
	mutate(lineage=case_when(lineage=="NV"~"Non-VOI/VOC",TRUE~lineage)) %>% 
	mutate(lineage=factor(lineage, c("Non-VOI/VOC","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>% 
	ggplot(aes(x=global_pars[["lod"]]-dp, col=lineage)) + 
		geom_density(size=1, alpha=0.6) + 
		theme_minimal() + 
		theme(text=element_text(size=14)) + 
		scale_color_manual(values=c("Non-VOI/VOC"="cornflowerblue","Alpha"="red2","Epsilon"="forestgreen","Delta"="orange","Unvaccinated"="black","Vaccinated"="darkblue")) + 
		labs(x="Ct",y="Density")


fig_delta_breakthrogh <- params_df_combined %>% 
	filter(lineage=="Delta") %>% 
	left_join(breakthroughdf, by="id") %>% 
	mutate(breakthrough=case_when(breakthrough=="No"~"Unvaccinated",breakthrough=="Yes"~"Breakthrough")) %>% 
	ggplot(aes(x=global_pars[["lod"]]-dp, fill=breakthrough)) + 
		geom_histogram(aes(y=..density..), binwidth=1, col="black", alpha=0.2, position="identity") +
		geom_density(alpha=0.2) + 
		facet_wrap(~breakthrough) + 
		theme_minimal() + 
		labs(title="Delta",x="Ct",y="Density") + 
		theme(text=element_text(size=14), legend.title=element_blank())

fig_minct_raw <- ct_dat_refined %>% 
	filter(!(B1351Status=="Yes" | P1Status=="Yes" | B1525Status=="Yes" | B1526Status=="Yes" | B15261Status=="Yes" | B15262Status=="Yes" | B16171Status=="Yes" | B16173Status=="Yes" | P2Status=="Yes")) %>% 
	group_by(PersonID) %>% 
	summarise(MinCt=min(CtT1), Lineage=first(Lineage), GreekLineage=first(GreekLineage)) %>% 
	mutate(GreekLineage=case_when(is.na(GreekLineage)~"Non-VOI/VOC",TRUE~GreekLineage)) %>% 
	mutate(GreekLineage=factor(GreekLineage, c("Non-VOI/VOC","Alpha","Epsilon","Delta"))) %>% 
	# mutate(Lineage=case_when(Lineage %in% c("B.1.1.7")~"Alpha",
	# Lineage%in%c("B.1.617.2","AY.1","AY.2","AY.3")~"Delta",
	# Lineage%in%c("B.1.427","B.1.429")~"Epsilon",
	# TRUE~"Other")) %>% 
	ggplot(aes(x=MinCt)) +
		geom_histogram(aes(y=..density..), binwidth=1, col="black", fill="lightgray") +
		# geom_density(adjust=1) + 
		theme_minimal() + 
		facet_wrap(~GreekLineage) + 
		labs(x="Minimum Ct", y="Density") + 
		theme(text=element_text(size=14)) 

# Look at raw traces in a fancy way: ------------------------------------------

thinsegmentsize <- 0.2
thinsegmentalpha <- 0.4
thicksegmentsize <- 1
thicksegmentalpha <- 0.8
smallpointsize <- 0.2
smallpointalpha <- 0.4
bigpointsize <- 0.75
bigpointalpha <- 1
nvcol <- "cornflowerblue"
alphacol <- "red2"
deltacol <- "orange"
unvaxcol <- "forestgreen"
vaxcol <- "darkblue"
pointcol <- "black"

ct_dat_refined_shifted <- ct_dat_refined %>% 
	select(PersonID, VaccineBreakthrough, GreekLineage, TestDateIndex, CtT1) %>%
	left_join((meanvalsindiv_combined %>% 
				group_by(id) %>% 
				slice(1) %>% 
				select(PersonID=id, tp)), by="PersonID") %>%  
	mutate(TestDateIndex=TestDateIndex-tp) %>% 
	select(-tp)

make_trajectory_plot_lin <- function(thislin, lincol="cornflowerblue", pointcol="black", thinsegmentsize=0.2, thinsegmentalpha=0.4, thicksegmentsize=1, thicksegmentalpha=0.8, smallpointsize=0.2, smallpointalpha=0.4, bigpointsize=0.75, bigpointalpha=1){

	highperson <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage == thislin) %>% 
		mutate(target=max(dp)) %>% 
		filter(dp==target) %>% 
		slice(1)

	lowperson <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage == thislin) %>% 
		mutate(target=min(dp)) %>% 
		filter(dp==target) %>% 
		slice(1)

	medperson <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage == thislin) %>% 
		mutate(target=median(dp)) %>% 
		mutate(targetdiff=abs(dp-target)) %>% 
		arrange(targetdiff) %>% 	
		slice(1)

	out <- meanvalsindiv_combined %>% 
	filter(lineage == thislin) %>% 
	ggplot() + 
		# Plot all the segments: 
		geom_segment(aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		geom_segment(aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		# Plot all the data: 
		geom_point(data=filter(ct_dat_refined_shifted, GreekLineage==thislin & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, size=smallpointsize, alpha=smallpointalpha) + 
		# Plot the high person lines: 
		geom_segment(data=highperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thicksegmentsize, alpha=1) + 
		geom_segment(data=highperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thicksegmentsize, alpha=1) + 
		# Plot the low person lines:
		geom_segment(data=lowperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thicksegmentsize, alpha=bigpointalpha) + 
		geom_segment(data=lowperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thicksegmentsize, alpha=bigpointalpha) + 
		# Plot the median person lines:
		geom_segment(data=medperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thicksegmentsize, alpha=bigpointalpha) + 
		geom_segment(data=medperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thicksegmentsize, alpha=bigpointalpha) + 
		# Plot high, low, and median points: 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID==highperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=24, alpha=bigpointalpha) + 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID==lowperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=25, alpha=bigpointalpha) + 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID==medperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=22, alpha=bigpointalpha) + 
		# Plot characteristics: 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
		theme_minimal() + 
		theme(text=element_text(size=18)) + 
		labs(x="Days since peak viral concentration", y="Ct") 

		return(out)

}



make_trajectory_plot_lin_bysample <- function(thislin, lincol="cornflowerblue", pointcol="black", thinsegmentsize=0.2, thinsegmentalpha=0.4, thicksegmentsize=1, thicksegmentalpha=0.8, smallpointsize=0.2, smallpointalpha=0.4, bigpointsize=0.75, bigpointalpha=1, ntrajectories=3){

	mostsampledids <- ct_dat_refined_shifted %>% 
		filter(GreekLineage==thislin) %>% 
		filter(CtT1<40) %>% 
		group_by(PersonID) %>% 
		summarise(NAbove=n()) %>% 
		arrange(desc(NAbove)) %>% 
		slice(1:ntrajectories) %>% 
		pull(PersonID) 

	mostsampled <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage==thislin) %>% 
		filter(id %in% mostsampledids)


	out <- meanvalsindiv_combined %>% 
	filter(lineage == thislin) %>% 
	ggplot() + 
		# Plot all the segments: 
		# geom_segment(aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		# geom_segment(aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		# Plot all the data: 
		geom_point(data=filter(ct_dat_refined_shifted, GreekLineage==thislin & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=lincol, size=smallpointsize, alpha=smallpointalpha) + 
		# Plot the special person trajectories: 
		geom_segment(data=mostsampled, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thicksegmentsize, alpha=1) + 
		geom_segment(data=mostsampled, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thicksegmentsize, alpha=1) + 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID%in%mostsampledids & CtT1<40), aes(x=TestDateIndex, y=CtT1, shape=factor(PersonID)), col=pointcol, fill=pointcol, size=bigpointsize, alpha=bigpointalpha) + 
		# Plot characteristics: 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
		theme_minimal() + 
		theme(text=element_text(size=18), legend.position="none") + 
		labs(x="Days since peak viral concentration", y="Ct") 

		return(out)

}




make_trajectory_plot_vax_bysample <- function(thisvax, lincol="cornflowerblue", pointcol="black", thinsegmentsize=0.2, thinsegmentalpha=0.4, thicksegmentsize=1, thicksegmentalpha=0.8, smallpointsize=0.2, smallpointalpha=0.4, bigpointsize=0.75, bigpointalpha=1, ntrajectories=3){

	thisvaxmarker <- ifelse(thisvax=="Vaccinated","Yes","No")

	mostsampledids <- ct_dat_refined_shifted %>% 
		filter(VaccineBreakthrough==thisvaxmarker) %>% 
		filter(CtT1<40) %>% 
		group_by(PersonID) %>% 
		summarise(NAbove=n()) %>% 
		arrange(desc(NAbove)) %>% 
		slice(1:ntrajectories) %>% 
		pull(PersonID) 

	mostsampled <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage==thisvax) %>% 
		filter(id %in% mostsampledids)


	out <- meanvalsindiv_combined %>% 
	filter(lineage == thisvax) %>% 
	ggplot() + 
		# Plot all the segments: 
		# geom_segment(aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		# geom_segment(aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		# Plot all the data: 
		geom_point(data=filter(ct_dat_refined_shifted, VaccineBreakthrough==thisvaxmarker & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=lincol, size=smallpointsize, alpha=smallpointalpha) + 
		# Plot the special person trajectories: 
		geom_segment(data=mostsampled, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=lincol, size=thicksegmentsize, alpha=1) + 
		geom_segment(data=mostsampled, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=lincol, size=thicksegmentsize, alpha=1) + 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID%in%mostsampledids & CtT1<40), aes(x=TestDateIndex, y=CtT1, shape=factor(PersonID)), col=pointcol, fill=pointcol, size=bigpointsize, alpha=bigpointalpha) + 
		# Plot characteristics: 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
		theme_minimal() + 
		theme(text=element_text(size=18), legend.position="none") + 
		labs(x="Days since peak viral concentration", y="Ct") 

		return(out)

}

make_raw_plot_vax_bysample <- function(thisvax, lincol="cornflowerblue", smallpointsize=0.2, smallpointalpha=0.4, tailpointalpha=0.4, xlims=c(-14,30), ylims=c(40,11)){

	thisvaxmarker <- ifelse(thisvax=="Vaccinated","Yes","No")

	out <- ct_dat_refined_shifted %>% 
		filter(VaccineBreakthrough==thisvaxmarker) %>% 
		filter(CtT1 < 40) %>% 
		left_join((meanvalsindiv_combined %>% 
					filter(lineage==thisvax) %>% 
					mutate(wr=wr-tp) %>% 
					select(PersonID=id, wr)),
				by="PersonID"
			) %>% 
		mutate(alpha=case_when(
			TestDateIndex <= wr ~ smallpointalpha,
			TRUE ~ tailpointalpha)) %>% 
		ggplot() + 
			# Plot all the data: 
			geom_point(aes(x=TestDateIndex, y=CtT1, alpha=alpha), col=lincol, size=smallpointsize) + 
			# Plot characteristics: 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml)), limits=ylims, breaks=c(40,35,30,25,20,15), labels=c("(-)","35","30","25","20","15")) + 
			theme_minimal() + 
			theme(text=element_text(size=18), legend.position="none") + 
			labs(x="Days since peak viral concentration", y="Ct") + 
			scale_x_continuous(limits=xlims) + 
			scale_alpha_identity() 

		return(out)

}

make_raw_plot_lin_bysample <- function(thislin, lincol="cornflowerblue", smallpointsize=0.2, smallpointalpha=0.4, tailpointalpha=0.4, xlims=c(-14,30), ylims=c(40,11)){

	out <- ct_dat_refined_shifted %>% 
		filter(GreekLineage==thislin) %>% 
		filter(CtT1 < 40) %>% 
		left_join((meanvalsindiv_combined %>% 
					filter(lineage==thislin) %>% 
					mutate(wr=wr-tp) %>% 
					select(PersonID=id, wr)),
				by="PersonID"
			) %>% 
		mutate(alpha=case_when(
			TestDateIndex <= wr ~ smallpointalpha,
			TRUE ~ tailpointalpha)) %>% 
		ggplot() + 
			# Plot all the data: 
			geom_point(aes(x=TestDateIndex, y=CtT1, alpha=alpha), col=lincol, size=smallpointsize) + 
			# Plot characteristics: 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml)), limits=ylims, breaks=c(40,35,30,25,20,15), labels=c("(-)","35","30","25","20","15")) + 
			theme_minimal() + 
			theme(text=element_text(size=18), legend.position="none") + 
			labs(x="Days since peak viral concentration", y="Ct") + 
			scale_x_continuous(limits=xlims) + 
			scale_alpha_identity() 

		return(out)

}

make_trajectory_plot_vax <- function(thisvax, lincol="cornflowerblue", pointcol="black", thinsegmentsize=0.2, thinsegmentalpha=0.4, thicksegmentsize=1, thicksegmentalpha=0.8, smallpointsize=0.2, smallpointalpha=0.4, bigpointsize=0.75, bigpointalpha=1){

	highperson <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage == thisvax) %>% 
		mutate(target=max(dp)) %>% 
		filter(dp==target) %>% 
		slice(1)

	lowperson <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage == thisvax) %>% 
		mutate(target=min(dp)) %>% 
		filter(dp==target) %>% 
		slice(1)

	medperson <- meanvalsindiv_combined %>% 
		ungroup() %>% 
		filter(lineage == thisvax) %>% 
		mutate(target=median(dp)) %>% 
		mutate(targetdiff=abs(dp-target)) %>% 
		arrange(targetdiff) %>% 	
		slice(1)

	thisvaxmarker <- ifelse(thisvax=="Vaccinated","Yes","No")

	out <- meanvalsindiv_combined %>% 
	filter(lineage == thisvax) %>% 
	ggplot() + 
		# Plot all the segments: 
		geom_segment(aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		geom_segment(aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
		# Plot all the data: 
		geom_point(data=filter(ct_dat_refined_shifted, VaccineBreakthrough==thisvaxmarker & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, size=smallpointsize, alpha=smallpointalpha) + 
		# Plot the high person lines: 
		geom_segment(data=highperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thicksegmentsize, alpha=1) + 
		geom_segment(data=highperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thicksegmentsize, alpha=1) + 
		# Plot the low person lines:
		geom_segment(data=lowperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
		geom_segment(data=lowperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
		# Plot the median person lines:
		geom_segment(data=medperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
		geom_segment(data=medperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
		# Plot high, low, and median points: 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID==highperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=24, alpha=bigpointalpha) + 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID==lowperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=25, alpha=bigpointalpha) + 
		geom_point(data=filter(ct_dat_refined_shifted, PersonID==medperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=22, alpha=bigpointalpha) + 
		# Plot characteristics: 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
		theme_minimal() + 
		theme(text=element_text(size=18)) + 
		labs(x="Days since peak viral concentration", y="Ct") 

		return(out)

}

fig_delta_sample_trajectories <- make_trajectory_plot_lin_bysample("Delta", lincol="orange", thicksegmentsize=0.5, thicksegmentalpha=0.5, bigpointsize=2, bigpointalpha=1, smallpointsize=0.5, smallpointalpha=1, ntrajectories=0)

fig_alpha_sample_trajectories <- make_trajectory_plot_lin_bysample("Alpha", lincol="red2", thicksegmentsize=0.5, thicksegmentalpha=0.5, bigpointsize=2, bigpointalpha=1, smallpointsize=0.5, smallpointalpha=1, ntrajectories=2)

fig_nv_sample_trajectories <- make_trajectory_plot_lin_bysample("NV", lincol="cornflowerblue", thicksegmentsize=0.5, thicksegmentalpha=0.5, bigpointsize=2, bigpointalpha=1, smallpointsize=0.5, smallpointalpha=1, ntrajectories=2)

fig_vax_sample_trajectories <- make_trajectory_plot_vax_bysample("Vaccinated", lincol="darkblue", thicksegmentsize=0.5, thicksegmentalpha=0.5, bigpointsize=2, bigpointalpha=1, smallpointsize=0.5, smallpointalpha=1, ntrajectories=2)

fig_unvax_sample_trajectories <- make_trajectory_plot_vax_bysample("Unvaccinated", lincol="forestgreen", thicksegmentsize=0.5, thicksegmentalpha=0.5, bigpointsize=2, bigpointalpha=1, smallpointsize=0.5, smallpointalpha=1, ntrajectories=2)




fig_delta_raw_points <- make_raw_plot_lin_bysample(thislin="Delta", lincol="orange", smallpointsize=1.5, smallpointalpha=1, tailpointalpha=0.3, xlims=c(-14,35))
fig_alpha_raw_points <- make_raw_plot_lin_bysample(thislin="Alpha", lincol="red2", smallpointsize=1.5, smallpointalpha=1, tailpointalpha=0.3, xlims=c(-14,35))
fig_nv_raw_points <- make_raw_plot_lin_bysample(thislin="NV", lincol="cornflowerblue", smallpointsize=1.5, smallpointalpha=1, tailpointalpha=0.3, xlims=c(-14,35))

fig_vax_raw_points <- make_raw_plot_vax_bysample(thisvax="Vaccinated", lincol="darkblue", smallpointsize=1.5, smallpointalpha=1, tailpointalpha=0.3, xlims=c(-14,35))
fig_unvax_raw_points <- make_raw_plot_vax_bysample(thisvax="Unvaccinated", lincol="forestgreen", smallpointsize=1.5, smallpointalpha=1, tailpointalpha=0.3, xlims=c(-14,35))




temp613 <- meanvalsindiv_combined %>% 	
	filter(lineage=="NV") %>% 
	filter(id==613) %>% 
	mutate(wp=wp-tp, wr=wr-tp) 

fig_example_alpha <- ct_dat_refined_shifted %>% 
	filter(PersonID==613) %>% 
	select(id=PersonID, t=TestDateIndex, y=CtT1) %>% 
	trim_negatives(global_pars) %>% 
	left_join((meanvalsindiv_combined %>% 
			filter(lineage=="NV") %>% 
			mutate(wr=wr-tp) %>% 
			select(id, wr)), 
	by="id") %>% 
	mutate(alpha=case_when(t <= wr ~ 1, TRUE ~ 0.3)) %>% 
	ggplot() + 
		geom_point(aes(x=t, y=y, alpha=alpha),size=4) + 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml)), breaks=seq(from=40, to=10, by=-5), minor_breaks=seq(from=40, to=10, by=-2.5), labels=c("(-)", "35","30","25","20","15","10")) + 
		geom_segment(data=temp613, aes(x=-Inf, xend=-wp, y=40, yend=40), size=1) + 
		geom_segment(data=temp613, aes(x=-wp, xend=-0, y=40, yend=40-dp), size=1) + 
		geom_segment(data=temp613, aes(x=0, xend=wr, y=40-dp, yend=40), size=1) + 
		geom_segment(data=temp613, aes(x=wr, xend=Inf, y=40, yend=40), size=1) + 
		theme_minimal() + 
		scale_alpha_identity() + 
		theme(text=element_text(size=24)) + 
		labs(x="Days since peak viral concentration", y="Ct")
	




# filter(ct_dat_refined, PersonID==medperson$id) %>% filter(CtT1<40) %>% select(PersonID, TestDateIndex, CtT1) 


# highperson <- meanvalsindiv_combined %>% 
# 	ungroup() %>% 
# 	filter(lineage == "Delta") %>% 
# 	mutate(target=max(dp)) %>% 
# 	filter(dp==target) %>% 
# 	slice(1)

# lowperson <- meanvalsindiv_combined %>% 
# 	ungroup() %>% 
# 	filter(lineage == "Delta") %>% 
# 	mutate(target=min(dp)) %>% 
# 	filter(dp==target) %>% 
# 	slice(1)

# medperson <- meanvalsindiv_combined %>% 
# 	ungroup() %>% 
# 	filter(lineage == "Delta") %>% 
# 	mutate(target=median(dp)) %>% 
# 	mutate(targetdiff=abs(dp-target)) %>% 
# 	arrange(targetdiff) %>% 	
# 	slice(1)

# meanvalsindiv_combined %>% 
# 	filter(lineage == "Delta") %>% 
# 	ggplot() + 
# 		# Plot all the segments: 
# 		geom_segment(aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
# 		geom_segment(aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thinsegmentsize, alpha=thinsegmentalpha) + 
# 		# Plot all the data: 
# 		geom_point(data=filter(ct_dat_refined_shifted, GreekLineage=="Delta" & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, size=smallpointsize, alpha=smallpointalpha) + 
# 		# Plot the high person lines: 
# 		geom_segment(data=highperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thicksegmentsize, alpha=1) + 
# 		geom_segment(data=highperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thicksegmentsize, alpha=1) + 
# 		# Plot the low person lines:
# 		geom_segment(data=lowperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
# 		geom_segment(data=lowperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
# 		# Plot the median person lines:
# 		geom_segment(data=medperson, aes(x=-wp, xend=0, y=global_pars[["lod"]], yend=global_pars[["lod"]]-dp), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
# 		geom_segment(data=medperson, aes(x=0, xend=wr, y=global_pars[["lod"]]-dp, yend=global_pars[["lod"]]), col=deltacol, size=thicksegmentsize, alpha=bigpointalpha) + 
# 		# Plot high, low, and median points: 
# 		geom_point(data=filter(ct_dat_refined_shifted, PersonID==highperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=24, alpha=bigpointalpha) + 
# 		geom_point(data=filter(ct_dat_refined_shifted, PersonID==lowperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=25, alpha=bigpointalpha) + 
# 		geom_point(data=filter(ct_dat_refined_shifted, PersonID==medperson$id & CtT1<40), aes(x=TestDateIndex, y=CtT1), col=pointcol, fill=pointcol, size=bigpointsize, shape=22, alpha=bigpointalpha) + 
# 		# Plot characteristics: 
# 		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
# 		theme_minimal() + 
# 		theme(text=element_text(size=14)) + 
# 		labs(x="Days since peak viral concentration", y="Ct") 

# Save: -----------------------------------------------------------------------
ggsave(fig_whiskers_dp, file=paste0(savedir_parent,"whiskers_dp.pdf"), width=6, height=4)
ggsave(fig_whiskers_wp, file=paste0(savedir_parent,"whiskers_wp.pdf"), width=6, height=4)
ggsave(fig_whiskers_wr, file=paste0(savedir_parent,"whiskers_wr.pdf"), width=6, height=4)
ggsave(fig_whiskers_infdur, file=paste0(savedir_parent,"whiskers_infdur.pdf"), width=6, height=4)

ggsave(fig_whiskers_dp_noepsilon, file=paste0(savedir_parent,"whiskers_dp_noepsilon.pdf"), width=6, height=4)
ggsave(fig_whiskers_wp_noepsilon, file=paste0(savedir_parent,"whiskers_wp_noepsilon.pdf"), width=6, height=4)
ggsave(fig_whiskers_wr_noepsilon, file=paste0(savedir_parent,"whiskers_wr_noepsilon.pdf"), width=6, height=4)
ggsave(fig_whiskers_infdur_noepsilon, file=paste0(savedir_parent,"whiskers_infdur_noepsilon.pdf"), width=6, height=4)

ggsave(fig_whiskers_dp_justvariants, file=paste0(savedir_parent,"whiskers_dp_justvariants.pdf"), width=6, height=4)
ggsave(fig_whiskers_wp_justvariants, file=paste0(savedir_parent,"whiskers_wp_justvariants.pdf"), width=6, height=4)
ggsave(fig_whiskers_wr_justvariants, file=paste0(savedir_parent,"whiskers_wr_justvariants.pdf"), width=6, height=4)
ggsave(fig_whiskers_infdur_justvariants, file=paste0(savedir_parent,"whiskers_infdur_justvariants.pdf"), width=6, height=4)

ggsave(fig_whiskers_dp_justvax, file=paste0(savedir_parent,"whiskers_dp_justvax.pdf"), width=6, height=4)
ggsave(fig_whiskers_wp_justvax, file=paste0(savedir_parent,"whiskers_wp_justvax.pdf"), width=6, height=4)
ggsave(fig_whiskers_wr_justvax, file=paste0(savedir_parent,"whiskers_wr_justvax.pdf"), width=6, height=4)
ggsave(fig_whiskers_infdur_justvax, file=paste0(savedir_parent,"whiskers_infdur_justvax.pdf"), width=6, height=4)

ggsave(fig_ct_trajectory_alpha, file=paste0(savedir_parent,"ct_trajectory_alpha.pdf"), width=8, height=5)
ggsave(fig_ct_trajectory_epsilon, file=paste0(savedir_parent,"ct_trajectory_epsilon.pdf"), width=8, height=5)
ggsave(fig_ct_trajectory_delta, file=paste0(savedir_parent,"ct_trajectory_delta.pdf"), width=8, height=5)
ggsave(fig_ct_trajectory_vax, file=paste0(savedir_parent,"ct_trajectory_vax.pdf"), width=8, height=5)

imap(figlist_ct_fit_nonVOIVOC, ~ggsave(.x,file=paste0(savedir_parent,"ct_fit_nonVOIVOC_",.y,".pdf"),width=8,height=8))
# ggsave(fig_ct_fit_nonVOIVOC, file=paste0(savedir_parent,"ct_fit_nonVOIVOC.pdf"), width=8, height=8)
ggsave(fig_ct_fit_alpha, file=paste0(savedir_parent,"ct_fit_alpha.pdf"), width=8, height=8)
ggsave(fig_ct_fit_epsilon, file=paste0(savedir_parent,"ct_fit_epsilon.pdf"), width=8, height=8)
ggsave(fig_ct_fit_delta, file=paste0(savedir_parent,"ct_fit_delta.pdf"), width=8, height=8)
imap(figlist_ct_fit_unvaccinated, ~ggsave(.x,file=paste0(savedir_parent,"ct_fit_unvaccinated_",.y,".pdf"),width=8,height=8))
imap(figlist_ct_fit_vaccinated, ~ggsave(.x,file=paste0(savedir_parent,"ct_fit_vaccinated_",.y,".pdf"),width=8,height=8))

ggsave(fig_lineage_hists, file=paste0(savedir_parent,"lineage_hists.pdf"), width=8, height=5)
ggsave(fig_lineage_hists_justvariants, file=paste0(savedir_parent,"lineage_hists_justvariants.pdf"), width=8, height=4)

ggsave(fig_delta_breakthrogh, file=paste0(savedir_parent,"delta_breakthrogh.pdf"), width=8, height=5)

ggsave(fig_minct_raw, file=paste0(savedir_parent,"minct_raw.pdf"), width=8, height=5)

ggsave(fig_delta_sample_trajectories, file=paste0(savedir_parent,"delta_sample_trajectories.pdf"), width=8, height=5)
ggsave(fig_alpha_sample_trajectories, file=paste0(savedir_parent,"alpha_sample_trajectories.pdf"), width=8, height=5)
ggsave(fig_nv_sample_trajectories, file=paste0(savedir_parent,"nv_sample_trajectories.pdf"), width=8, height=5)
ggsave(fig_vax_sample_trajectories, file=paste0(savedir_parent,"vax_sample_trajectories.pdf"), width=8, height=5)
ggsave(fig_unvax_sample_trajectories, file=paste0(savedir_parent,"unvax_sample_trajectories.pdf"), width=8, height=5)

ggsave(fig_delta_raw_points, file=paste0(savedir_parent,"delta_raw_points.pdf"), width=8, height=5)
ggsave(fig_alpha_raw_points, file=paste0(savedir_parent,"alpha_raw_points.pdf"), width=8, height=5)
ggsave(fig_nv_raw_points, file=paste0(savedir_parent,"nv_raw_points.pdf"), width=8, height=5)
ggsave(fig_vax_raw_points, file=paste0(savedir_parent,"vax_raw_points.pdf"), width=8, height=5)
ggsave(fig_unvax_raw_points, file=paste0(savedir_parent,"unvax_raw_points.pdf"), width=8, height=5)
ggsave(fig_example_alpha, file=paste0(savedir_parent,"example_alpha.pdf"), width=8, height=5)

# Save png: --------------------------------------------------------------------
ggsave(fig_whiskers_dp, file=paste0(savedir_parent,"whiskers_dp.png"), width=6, height=4)
ggsave(fig_whiskers_wp, file=paste0(savedir_parent,"whiskers_wp.png"), width=6, height=4)
ggsave(fig_whiskers_wr, file=paste0(savedir_parent,"whiskers_wr.png"), width=6, height=4)
ggsave(fig_whiskers_infdur, file=paste0(savedir_parent,"whiskers_infdur.png"), width=6, height=4)

ggsave(fig_whiskers_dp_noepsilon, file=paste0(savedir_parent,"whiskers_dp_noepsilon.png"), width=6, height=4)
ggsave(fig_whiskers_wp_noepsilon, file=paste0(savedir_parent,"whiskers_wp_noepsilon.png"), width=6, height=4)
ggsave(fig_whiskers_wr_noepsilon, file=paste0(savedir_parent,"whiskers_wr_noepsilon.png"), width=6, height=4)
ggsave(fig_whiskers_infdur_noepsilon, file=paste0(savedir_parent,"whiskers_infdur_noepsilon.png"), width=6, height=4)

ggsave(fig_whiskers_dp_justvariants, file=paste0(savedir_parent,"whiskers_dp_justvariants.png"), width=6, height=4)
ggsave(fig_whiskers_wp_justvariants, file=paste0(savedir_parent,"whiskers_wp_justvariants.png"), width=6, height=4)
ggsave(fig_whiskers_wr_justvariants, file=paste0(savedir_parent,"whiskers_wr_justvariants.png"), width=6, height=4)
ggsave(fig_whiskers_infdur_justvariants, file=paste0(savedir_parent,"whiskers_infdur_justvariants.png"), width=6, height=4)

ggsave(fig_whiskers_dp_justvax, file=paste0(savedir_parent,"whiskers_dp_justvax.png"), width=6, height=4)
ggsave(fig_whiskers_wp_justvax, file=paste0(savedir_parent,"whiskers_wp_justvax.png"), width=6, height=4)
ggsave(fig_whiskers_wr_justvax, file=paste0(savedir_parent,"whiskers_wr_justvax.png"), width=6, height=4)
ggsave(fig_whiskers_infdur_justvax, file=paste0(savedir_parent,"whiskers_infdur_justvax.png"), width=6, height=4)

ggsave(fig_ct_trajectory_alpha, file=paste0(savedir_parent,"ct_trajectory_alpha.png"), width=8, height=5)
ggsave(fig_ct_trajectory_epsilon, file=paste0(savedir_parent,"ct_trajectory_epsilon.png"), width=8, height=5)
ggsave(fig_ct_trajectory_delta, file=paste0(savedir_parent,"ct_trajectory_delta.png"), width=8, height=5)
ggsave(fig_ct_trajectory_vax, file=paste0(savedir_parent,"ct_trajectory_vax.png"), width=8, height=5)

imap(figlist_ct_fit_nonVOIVOC, ~ggsave(.x,file=paste0(savedir_parent,"ct_fit_nonVOIVOC_",.y,".png"),width=8,height=8))
# ggsave(fig_ct_fit_nonVOIVOC, file=paste0(savedir_parent,"ct_fit_nonVOIVOC.png"), width=8, height=8)
ggsave(fig_ct_fit_alpha, file=paste0(savedir_parent,"ct_fit_alpha.png"), width=8, height=8)
ggsave(fig_ct_fit_epsilon, file=paste0(savedir_parent,"ct_fit_epsilon.png"), width=8, height=8)
ggsave(fig_ct_fit_delta, file=paste0(savedir_parent,"ct_fit_delta.png"), width=8, height=8)
imap(figlist_ct_fit_unvaccinated, ~ggsave(.x,file=paste0(savedir_parent,"ct_fit_unvaccinated_",.y,".png"),width=8,height=8))
imap(figlist_ct_fit_vaccinated, ~ggsave(.x,file=paste0(savedir_parent,"ct_fit_vaccinated_",.y,".png"),width=8,height=8))

ggsave(fig_lineage_hists, file=paste0(savedir_parent,"lineage_hists.png"), width=8, height=5)
ggsave(fig_lineage_hists_justvariants, file=paste0(savedir_parent,"lineage_hists_justvariants.png"), width=8, height=4)

ggsave(fig_delta_breakthrogh, file=paste0(savedir_parent,"delta_breakthrogh.png"), width=8, height=5)

ggsave(fig_minct_raw, file=paste0(savedir_parent,"minct_raw.png"), width=8, height=5)

ggsave(fig_delta_sample_trajectories, file=paste0(savedir_parent,"delta_sample_trajectories.png"), width=8, height=5)
ggsave(fig_alpha_sample_trajectories, file=paste0(savedir_parent,"alpha_sample_trajectories.png"), width=8, height=5)
ggsave(fig_nv_sample_trajectories, file=paste0(savedir_parent,"nv_sample_trajectories.png"), width=8, height=5)
ggsave(fig_vax_sample_trajectories, file=paste0(savedir_parent,"vax_sample_trajectories.png"), width=8, height=5)
ggsave(fig_unvax_sample_trajectories, file=paste0(savedir_parent,"unvax_sample_trajectories.png"), width=8, height=5)

ggsave(fig_delta_raw_points, file=paste0(savedir_parent,"delta_raw_points.png"), width=8, height=5)
ggsave(fig_alpha_raw_points, file=paste0(savedir_parent,"alpha_raw_points.png"), width=8, height=5)
ggsave(fig_nv_raw_points, file=paste0(savedir_parent,"nv_raw_points.png"), width=8, height=5)
ggsave(fig_vax_raw_points, file=paste0(savedir_parent,"vax_raw_points.png"), width=8, height=5)
ggsave(fig_unvax_raw_points, file=paste0(savedir_parent,"unvax_raw_points.png"), width=8, height=5)
ggsave(fig_example_alpha, file=paste0(savedir_parent,"example_alpha.png"), width=8, height=5)

# Some KS tests: --------------------------------------------------------------
# Alpha:
ks_alpha_dp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Alpha")$dp,
	filter(meanvalsindiv_combined,lineage=="NV")$dp
	)$p.value

ks_alpha_wp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Alpha")$wp,
	filter(meanvalsindiv_combined,lineage=="NV")$wp
	)$p.value

ks_alpha_wr <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Alpha")$wr,
	filter(meanvalsindiv_combined,lineage=="NV")$wr
	)$p.value

ks_alpha_infdur <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Alpha")$infdur,
	filter(meanvalsindiv_combined,lineage=="NV")$infdur
	)$p.value

c(ks_alpha_dp, ks_alpha_wp, ks_alpha_wr, ks_alpha_infdur)

# Epsilon:
ks_epsilon_dp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Epsilon")$dp,
	filter(meanvalsindiv_combined,lineage=="NV")$dp
	)$p.value

ks_epsilon_wp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Epsilon")$wp,
	filter(meanvalsindiv_combined,lineage=="NV")$wp
	)$p.value

ks_epsilon_wr <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Epsilon")$wr,
	filter(meanvalsindiv_combined,lineage=="NV")$wr
	)$p.value

ks_epsilon_infdur <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Epsilon")$infdur,
	filter(meanvalsindiv_combined,lineage=="NV")$infdur
	)$p.value

c(ks_epsilon_dp, ks_epsilon_wp, ks_epsilon_wr, ks_epsilon_infdur)

# Delta:
ks_delta_dp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Delta")$dp,
	filter(meanvalsindiv_combined,lineage=="NV")$dp
	)$p.value

ks_delta_wp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Delta")$wp,
	filter(meanvalsindiv_combined,lineage=="NV")$wp
	)$p.value

ks_delta_wr <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Delta")$wr,
	filter(meanvalsindiv_combined,lineage=="NV")$wr
	)$p.value

ks_delta_infdur <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Delta")$infdur,
	filter(meanvalsindiv_combined,lineage=="NV")$infdur
	)$p.value

c(ks_delta_dp, ks_delta_wp, ks_delta_wr, ks_delta_infdur)

# Vaccination:
ks_vax_dp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Vaccinated")$dp,
	filter(meanvalsindiv_combined,lineage=="Unvaccinated")$dp
	)$p.value

ks_vax_wp <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Vaccinated")$wp,
	filter(meanvalsindiv_combined,lineage=="Unvaccinated")$wp
	)$p.value

ks_vax_wr <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Vaccinated")$wr,
	filter(meanvalsindiv_combined,lineage=="Unvaccinated")$wr
	)$p.value

ks_vax_infdur <- ks.test(
	filter(meanvalsindiv_combined,lineage=="Vaccinated")$infdur,
	filter(meanvalsindiv_combined,lineage=="Unvaccinated")$infdur
	)$p.value

c(ks_vax_dp, ks_vax_wp, ks_vax_wr, ks_vax_infdur)


ks_df <- tibble(
	Comparison=c("AlphaVsNV","EpsilonVsNV","DeltaVsNV","VaxVsUnvax"),
	Peak=c(ks_alpha_dp, ks_epsilon_dp, ks_delta_dp, ks_vax_dp),
	Proliferation=c(ks_alpha_wp, ks_epsilon_wp, ks_delta_wp, ks_vax_wp),
	Clearance=c(ks_alpha_wr, ks_epsilon_wr, ks_delta_wr, ks_vax_wr),
	Duration=c(ks_alpha_infdur, ks_epsilon_infdur, ks_delta_infdur, ks_vax_infdur),
	) %>% 
	mutate(Peak=round(Peak,3)) %>% 
	mutate(Proliferation=round(Proliferation,3)) %>% 
	mutate(Clearance=round(Clearance,3)) %>% 
	mutate(Duration=round(Duration,3))

# Summarize individual values: ------------------------------------------------

indiv_table_ct <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		min=min(global_pars[["lod"]]-dp),
		lwr=quantile(global_pars[["lod"]]-dp, 0.25),
		median=median(global_pars[["lod"]]-dp),
		upr=quantile(global_pars[["lod"]]-dp, 0.75),
		max=max(global_pars[["lod"]]-dp)
		) %>% 
	mutate(lineage=factor(lineage, c("NV","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>%
	arrange(lineage) %>%
	mutate(median=as.character(round(median,1))) %>%
	mutate(lwr=paste0(" [",as.character(round(lwr,1)),",")) %>%
	mutate(upr=paste0(" ",as.character(round(upr,1)),"]")) %>%
	mutate(min=paste0(" [",as.character(round(min,1)),",")) %>%
	mutate(max=paste0(" ",as.character(round(max,1)),"]")) %>%
	unite(summary, median, lwr, upr, min, max, sep="")

indiv_table_geml <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		min=min(convert_Ct_logGEML(global_pars[["lod"]]-dp)),
		lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dp), 0.25),
		median=median(convert_Ct_logGEML(global_pars[["lod"]]-dp)),
		upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dp), 0.75),
		max=max(convert_Ct_logGEML(global_pars[["lod"]]-dp))
		) %>% 
	mutate(lineage=factor(lineage, c("NV","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>%
	arrange(lineage) %>%
	mutate(median=as.character(round(median,1))) %>%
	mutate(lwr=paste0(" [",as.character(round(lwr,1)),",")) %>%
	mutate(upr=paste0(" ",as.character(round(upr,1)),"]")) %>%
	mutate(min=paste0(" [",as.character(round(min,1)),",")) %>%
	mutate(max=paste0(" ",as.character(round(max,1)),"]")) %>%
	unite(summary, median, lwr, upr, min, max, sep="")

indiv_table_wp <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		min=min(wp),
		lwr=quantile(wp, 0.25),
		median=median(wp),
		upr=quantile(wp, 0.75),
		max=max(wp)
		) %>% 
	mutate(lineage=factor(lineage, c("NV","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>%
	arrange(lineage) %>%
	mutate(median=as.character(round(median,1))) %>%
	mutate(lwr=paste0(" [",as.character(round(lwr,1)),",")) %>%
	mutate(upr=paste0(" ",as.character(round(upr,1)),"]")) %>%
	mutate(min=paste0(" [",as.character(round(min,1)),",")) %>%
	mutate(max=paste0(" ",as.character(round(max,1)),"]")) %>%
	unite(summary, median, lwr, upr, min, max, sep="")

indiv_table_wr <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		min=min(wr),
		lwr=quantile(wr, 0.25),
		median=median(wr),
		upr=quantile(wr, 0.75),
		max=max(wr)
		) %>% 
	mutate(lineage=factor(lineage, c("NV","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>%
	arrange(lineage) %>%
	mutate(median=as.character(round(median,1))) %>%
	mutate(lwr=paste0(" [",as.character(round(lwr,1)),",")) %>%
	mutate(upr=paste0(" ",as.character(round(upr,1)),"]")) %>%
	mutate(min=paste0(" [",as.character(round(min,1)),",")) %>%
	mutate(max=paste0(" ",as.character(round(max,1)),"]")) %>%
	unite(summary, median, lwr, upr, min, max, sep="")

indiv_table_infdur <- meanvalsindiv_combined %>% 
	group_by(lineage) %>% 
	summarise(
		min=min(infdur),
		lwr=quantile(infdur, 0.25),
		median=median(infdur),
		upr=quantile(infdur, 0.75),
		max=max(infdur)
		) %>% 
	mutate(lineage=factor(lineage, c("NV","Alpha","Epsilon","Delta","Unvaccinated","Vaccinated"))) %>%
	arrange(lineage) %>%
	mutate(median=as.character(round(median,1))) %>%
	mutate(lwr=paste0(" [",as.character(round(lwr,1)),",")) %>%
	mutate(upr=paste0(" ",as.character(round(upr,1)),"]")) %>%
	mutate(min=paste0(" [",as.character(round(min,1)),",")) %>%
	mutate(max=paste0(" ",as.character(round(max,1)),"]")) %>%
	unite(summary, median, lwr, upr, min, max, sep="")

print(indiv_table_ct)
print(indiv_table_geml)
print(indiv_table_wp)
print(indiv_table_wr)
print(indiv_table_infdur)


# Summarise population values: ------------------------------


format_dp <- function(df){
	addn_Ct <- df %>% 
		filter(statistic=="dp") %>% 
		mutate(statistic="MinCt") %>% 
		mutate(mean=global_pars[["lod"]]-mean) %>% 
		mutate(lwr=global_pars[["lod"]]-upr) %>% 
		mutate(upr=global_pars[["lod"]]-lwr)

	addn_geml <- addn_Ct %>% 
		mutate(statistic="MinGeml") %>% 
		mutate(mean=convert_Ct_logGEML(mean)) %>% 
		mutate(lwr=convert_Ct_logGEML(upr)) %>% 
		mutate(upr=convert_Ct_logGEML(lwr))

	out <- df %>% 
		filter(statistic!="dp") %>% 
		bind_rows(addn_Ct) %>% 
		bind_rows(addn_geml)
}

popmeantab <- shared_params_df_summary_forplotting %>% 
	format_dp() %>% 
	mutate(mean=as.character(round(mean,1))) %>% 
	mutate(lwr=paste0(" [",as.character(round(lwr,1)),",")) %>% 
	mutate(upr=paste0(" ",as.character(round(upr,1)),"]")) %>% 
	unite(interval, mean, lwr, upr, sep="") %>%
	select(-xval) %>% 
	pivot_wider(names_from=statistic, values_from=interval) %>%
	arrange(lineage) %>% 
	select(lineage, MinCt, MinGeml, wp, wr, id)
write_csv(popmeantab, file=paste0(savedir_parent,"popmeantab.csv"))


ct_dat_refined %>% 
	filter(CtT1 < 40) %>% 
	group_by(PersonID) %>% 
	summarise(NSamp=n()) %>% 
	summarise(min(NSamp), quantile(NSamp, 0.25), quantile(NSamp, 0.5), quantile(NSamp, 0.75), max(NSamp))

# system("say Done with saving figures overall")