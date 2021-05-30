# =============================================================================
# Fit model for all parameter sets
# =============================================================================

library(tidyverse) 
library(scales)

source('code/utils.R')
# NOTE: utils_private.R contains just one line specifying a directory on an external drive for storing the MCMC fits, since these can get large. The line looks like: 
# extdrive <- "/External/Drive/Path/"
# New users of this code should define their own utils_private.R script with a custom path. 
source('code/utils_private.R')
source("code/set_global_pars.R")
ct_dat_refined <- read_csv("data/ct_dat_refined.csv")
source("code/set_run_pars.R")

for(run_pars_index in c(1,3)){ # 1:length(run_pars_list)

	run_pars <- run_pars_list[[run_pars_index]]

	source("code/fit_posteriors_preamble.R")
	source("code/fit_posteriors.R")
	source("code/extract_parameters.R")

	source("code/make_figures.R")
	source("code/summarise_dists.R")

	# Set directory for saving figures: 
	savedir <- paste0("figures/run_pars_",run_pars_index,"/")
	# Set directory for saving MCMC output: 
	savedir_output <- paste0("output/run_pars_",run_pars_index,"/")

	source("code/save_figures.R")
	write_csv(dist_summary, file=paste0(savedir,"dist_summary.csv"))
	save(ct_fit, file=paste0(extdrive,savedir_output,"ct_fit.RData"))

	source('code/diagnose_fit.R')
	write_csv(warningtab, file=paste0(savedir,"warningtab.csv"))

	print(paste0("Finished parameter set ",run_pars_index,": ",Sys.time()))

}

source("code/make_figures_overall.R")

# =============================================================================
# Post-hoc figure making
# =============================================================================

# for(run_pars_index in 1:length(run_pars_list)){ 

# 	run_pars <- run_pars_list[[run_pars_index]]

# 	savedir <- paste0("figures/run_pars_",run_pars_index,"/")
# 	savedir_output <- paste0("output/run_pars_",run_pars_index,"/")

# 	source('code/fit_posteriors_preamble.R')
# 	load(file=paste0(extdrive,savedir_output,"ct_fit.RData"))
# 	source('code/extract_parameters.R')

# 	source("code/make_figures.R")
# 	source("code/summarise_dists.R")
# 	source("code/save_figures.R")
# 	write_csv(dist_summary, file=paste0(savedir,"dist_summary.csv"))

# 	source('code/diagnose_fit.R')
# 	write_csv(warningtab, file=paste0(savedir,"warningtab.csv"))

# 	print(paste0("Done with iteration ",run_pars_index))

# }

# launch_shinystan_nonblocking(ct_fit)
