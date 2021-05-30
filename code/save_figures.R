
ggsave(fig_ct_fit, file=paste0(savedir,"ct_fit.pdf"), width=8, height=8)
ggsave(fig_ct_fit, file=paste0(savedir,"ct_fit.png"), width=8, height=8)
ggsave(fig_dpmean, file=paste0(savedir,"dp_mean.pdf"), width=8, height=5)
ggsave(fig_dpmean_withprior, file=paste0(savedir,"dp_mean_withprior.pdf"), width=8, height=5)
ggsave(fig_gemlmean, file=paste0(savedir,"geml_mean.pdf"), width=8, height=5)
ggsave(fig_wpmean, file=paste0(savedir,"wp_mean.pdf"), width=8, height=5)
ggsave(fig_wpmean_withprior, file=paste0(savedir,"wp_mean_withprior.pdf"), width=8, height=5)
if(run_pars$trapfit==TRUE){
	ggsave(fig_wxmean, file=paste0(savedir,"wx_mean.pdf"), width=8, height=5)
	ggsave(fig_wxmean_withprior, file=paste0(savedir,"wx_mean_withprior.pdf"), width=8, height=5)
}
ggsave(fig_wrmean, file=paste0(savedir,"wr_mean.pdf"), width=8, height=5)
ggsave(fig_wrmean_withprior, file=paste0(savedir,"wr_mean_withprior.pdf"), width=8, height=5)
ggsave(fig_infdurmean, file=paste0(savedir,"infdur_mean.pdf"), width=8, height=5)
ggsave(fig_ct_trajectory_inference, file=paste0(savedir,"ct_trajectory_inference.pdf"), width=8, height=5)
ggsave(fig_indiv_trajectories, file=paste0(savedir,"indiv_trajectories.pdf"), width=8, height=5)
ggsave(fig_indiv_trajectories, file=paste0(savedir,"indiv_trajectories.png"), width=8, height=5)
ggsave(fig_indiv_trajectories_green, file=paste0(savedir,"indiv_trajectories_green.png"), width=8, height=5)
