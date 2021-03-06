library(tidyverse)
library(lazyeval)
library(rstan) 
library(shinystan) 
library(purrr)
library(data.table)
options(mc.cores=parallel::detectCores())
source('code/utils.R')
source("code/set_global_pars.R")


# Define a pared-down dataset for passing to Stan: 
# indiv_data <- ct_dat_refined %>% 
# 	filter(!(RowID %in% run_pars$excluded_rows)) %>%
# 	mutate(Special=case_when(RowID %in% run_pars$special_rows~1, TRUE~0)) %>%
# 	clean_person_id() %>%
# 	select(PersonID, PersonIDClean, TestDateIndex, CtT1, Special) %>%
# 	rename(id=PersonID) %>%
# 	rename(id_clean=PersonIDClean) %>% 
# 	rename(t=TestDateIndex) %>%
# 	rename(y=CtT1) %>%
# 	rename(special=Special) %>%
# 	trim_negatives(global_pars)
indiv_data <- ct_dat_refined %>% 
	filter(!(RowID %in% run_pars$excluded_rows)) %>%
	mutate(Special=case_when(RowID %in% run_pars$special_rows~1, TRUE~0)) %>%
	mutate(id=PersonID) %>%
	mutate(t=TestDateIndex) %>%
	mutate(y=CtT1) %>%
	mutate(special=Special) %>%
	trim_negatives(global_pars) %>% 
	clean_person_id() %>%
	mutate(id_clean=PersonIDClean) %>% 
	select(id, id_clean, t, y, special)

# Store the number of people we've kept: 
n_indiv <- length(unique(indiv_data$id))

special <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, special) %>%
	arrange(id) %>%
	pull(special)

# Useful dataframe for mapping official ID to Stan ID:
id_map <- indiv_data %>% 
	group_by(id) %>%
	summarise(id_clean=first(id_clean)) %>% 
	select(id, id_clean) %>%
	mutate(id_clean=as.character(id_clean))

# Useful dataframe for mapping PersonID to the 'special' entries
special_map <- indiv_data %>% 
	group_by(id) %>%
	summarise(special=first(special))

prior_pars <- list(
	special=special,
	tpsd=run_pars$tpsd,
	dpmin=run_pars$dpmin,
	dpmean_prior=run_pars$dpmean_prior,
	dpsd_prior=run_pars$dpsd_prior,
	wpmin=run_pars$wpmin,
	wpmax=run_pars$wpmax,
	wpmean_prior=run_pars$wpmean_prior,
	wpsd_prior=run_pars$wpsd_prior,
	wrmin=run_pars$wrmin,
	wrmax=run_pars$wrmax,
	wrmean_prior=run_pars$wrmean_prior,
	wrsd_prior=run_pars$wrsd_prior,
	sigma_max=run_pars$sigma_max,
	sigma_prior_scale=run_pars$sigma_prior_scale,
	lambda=run_pars$lambda,
	fpmean=run_pars$fpmean		# so that 90% of mass is <1 and 99% is <2
)	

if(run_pars$trapfit==1){
	prior_pars$wxmin <- run_pars$wxmin
	prior_pars$wxmax <- run_pars$wxmax
	prior_pars$wxmean_prior <- run_pars$wxmean_prior
	prior_pars$wxsd_prior <- run_pars$wxsd_prior
}
