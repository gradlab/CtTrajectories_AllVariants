library(tidyverse) 

# Re-assign Person IDs as a running index of 1:npeople
clean_person_id <- function(ct_dat){
	id_list <- sort(unique(ct_dat$PersonID))
	id_df <- data.frame(PersonID=id_list, PersonIDClean=1:length(id_list), stringsAsFactors=FALSE)
	out <- left_join(ct_dat, id_df, by="PersonID")
	return(out)
}

# Script for saving ggplot output stored in lists
save_figlist <- function(obj,dir,name,driver,width=8,height=5){
	mapply(function(x,y) ggsave(x,
		file=paste0(dir,name,"_",y,".",driver), width=width, height=height),
	obj,1:length(obj))
}

# Converts raw test dates to elapsed days since peak Ct:
make_test_date_index <- function(ct_dat){
	minctdf <- ct_dat %>% 
		arrange(PersonID, TestDate) %>%
		group_by(PersonID) %>%
		mutate(index=1:n()) %>%
		filter(CtT1==min(CtT1)) %>%
		slice(1) %>%
		select(PersonID, TestDate) %>%
		rename(MinCtTestDate=TestDate)

	out <- ct_dat %>% 
		left_join(minctdf, by="PersonID") %>%
		mutate(TestDateIndex=as.numeric(difftime(TestDate,MinCtTestDate,units="days"))) %>%
		select(-MinCtTestDate)

	return(out)
}

# Trims consecutive sequences of 3 or more test results at the limit of detection: 
trim_negatives <- function(indiv_data, global_pars){
	# lod <- 40
	out <- indiv_data %>% 
		split(.$id) %>% 
		map(~ arrange(., t)) %>% 
		map(~ mutate(., rowindex=1:n())) %>% 
		map(~ mutate(., ispositive=case_when(y<global_pars[["lod"]] ~ 1, TRUE~0))) %>% 
		map(~ mutate(., ispositive_lag=lag(ispositive))) %>%
		map(~ mutate(., ispositive_lag2=lag(ispositive,2))) %>%
		map(~ mutate(., ispositive_lead=lead(ispositive))) %>%
		map(~ mutate(., ispositive_lead2=lead(ispositive,2))) %>%
		map(~ filter(., ispositive==1 | ispositive_lag==1 | ispositive_lag2==1 | ispositive_lead==1 | ispositive_lead2==1)) %>%
		map(~ select(., id, id_clean, t, y, special)) %>%
		bind_rows()
	return(out)
}

# Laaunch shinystan without blocking the console: 
launch_shinystan_nonblocking <- function(fit) {
  library(future)
  plan(multisession)
  future(
    launch_shinystan(fit) #You can replace this with any other Shiny app
  )
}

# Like dplyr's sample_n, but sampling by group instead of by row: 
sample_n_groups <- function(grouped_df, size, replace=FALSE, weight=NULL){
	# From https://cmdlinetips.com/2019/07/how-to-randomly-select-groups-in-r-with-dplyr/
	grp_var <- grouped_df %>% 
	    groups %>%
	    unlist %>% 
	    as.character

	random_grp <- grouped_df %>% 
	    summarise() %>% 
	    slice_sample(n=size, replace=replace, weight_by=weight)# %>% 
	    # mutate(unique_id = 1:n())

	out <- grouped_df %>% 
	    right_join(random_grp, by=grp_var)

	return(out)
}

# For generating names for a matrix-turned-data frame: 
makenames <- function(parname, n_indiv){
	unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}

# For parsing the 'parameters' output from Stan: 
parseparam <- function(extracted_params, parname, n_indiv){
	as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}

# For collecting MCMC output into a useful data frame
make_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value") %>%
		select(-iteration) 
}

# For collecting MCMC output into a useful data frame (individual-level parameters)
make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value")
}

# For collecting MCMC output into a useful data frame (shared parameters)
make_shared_params_df <- function(extracted_params, parnames){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) 
		as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
		), cbind) %>%
		as_tibble() %>%
		mutate(iteration=1:n())
}

# Plot individual-level fitted Ct trajectories for everyone, with raw data: 
plot_ct_fit <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=100){
	with(as.list(global_pars),{
	params_df %>% 
		mutate(id_clean=as.integer(id_clean)) %>%
		group_by(id_clean) %>%
		sample_n(ntraces) %>% 
		ungroup() %>%
		ggplot() + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=indiv_data, aes(x=t, y=y,col=factor(special)), size=0.5) + 
			scale_color_manual(values=c("1"="red","0"="blue"), labels=c("1"="Variant","0"="Non-variant")) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
			labs(x="Time since min Ct (days)", y="Ct") + 
			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			facet_wrap(~id) # Change to id_clean to obscure identities
			})
}

# Turn off the background grid:
grid_off <- list(theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))

# Turn off y ticks:
y_ticks_off <- list(theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()))

# Convert Ct values to RNA copies per ml:
convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
	out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
	return(out) 
}

# Function used for plotting mean proliferation envelope:
srise <- function(x, dp, wp){
	out <- dp*(1+x/wp)
	return(out)
}
# Function used for plotting mean clearance envelope:
sfall <- function(x, dp, wr){
	out <- dp*(1-x/wr)
	return(out)
}

make_sample_trajectory_special <- function(shared_params_df, global_pars, siglevel=0.9, ge=FALSE, referencefill="blue", specialfill="red"){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wrmeanW, -wrmeanB) %>%
		rename(dpmeanW=dpmeanW_trans,
		dpmeanB=dpmeanB_trans,
		wpmeanW=wpmeanW_trans,
		wpmeanB=wpmeanB_trans,
		wrmeanW=wrmeanW_trans,
		wrmeanB=wrmeanB_trans)

	wp_mean_W <- mean(shared_params_df$wpmeanW)
	wp_lwr_W <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
	wp_upr_W <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

	wr_mean_W <- mean(shared_params_df$wrmeanW)
	wr_lwr_W <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
	wr_upr_W <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

	dp_mean_W <- mean(shared_params_df$dpmeanW)
	dp_lwr_W <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
	dp_upr_W <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

	xvals_proliferation_W <- seq(from=-wp_upr_W, 0, length.out=500)
	xvals_clearance_W <- seq(from=0, wr_upr_W, length.out=500)

	yvals_upr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
	yvals_upr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
	yvals_lwr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))

	# For symptomatic:
	wp_mean_B <- mean(shared_params_df$wpmeanB)
	wp_lwr_B <- quantile(shared_params_df$wpmeanB,(1-siglevel)/2)
	wp_upr_B <- quantile(shared_params_df$wpmeanB,1-(1-siglevel)/2)

	wr_mean_B <- mean(shared_params_df$wrmeanB)
	wr_lwr_B <- quantile(shared_params_df$wrmeanB,(1-siglevel)/2)
	wr_upr_B <- quantile(shared_params_df$wrmeanB,1-(1-siglevel)/2)

	dp_mean_B <- mean(shared_params_df$dpmeanB)
	dp_lwr_B <- quantile(shared_params_df$dpmeanB,(1-siglevel)/2)
	dp_upr_B <- quantile(shared_params_df$dpmeanB,1-(1-siglevel)/2)

	xvals_proliferation_B <- seq(from=-wp_upr_B, 0, length.out=500)
	xvals_clearance_B <- seq(from=0, wr_upr_B, length.out=500)

	yvals_upr_proliferation_B <- unlist(lapply(xvals_proliferation_B, 
	function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_B <- unlist(lapply(xvals_proliferation_B, 
	function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),(1-siglevel)/2)))
	yvals_upr_clearance_B <- unlist(lapply(xvals_clearance_B, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),1-(1-siglevel)/2)))
	yvals_lwr_clearance_B <- unlist(lapply(xvals_clearance_B, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),(1-siglevel)/2)))

	if(ge==FALSE){
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=yvals_lwr_proliferation_W, 
					yvals_upr=yvals_upr_proliferation_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=lod,yend=lod-dp_mean_W),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_B, 
					yvals_lwr=yvals_lwr_proliferation_B, 
					yvals_upr=yvals_upr_proliferation_B),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=-wp_mean_B,xend=0,y=lod,yend=lod-dp_mean_B),col=specialfill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=yvals_lwr_clearance_W, 
					yvals_upr=yvals_upr_clearance_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=lod-dp_mean_W,yend=lod),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_B, 
					yvals_lwr=yvals_lwr_clearance_B, 
					yvals_upr=yvals_upr_clearance_B),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=0,xend=wr_mean_B,y=lod-dp_mean_B,yend=lod),col=specialfill) + 
			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Ct") + 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + # "log"["10"]" RNA copies/ml"
			theme(text=element_text(size=18))
	} else {
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=(yvals_lwr_proliferation_W), 
					yvals_upr=(yvals_upr_proliferation_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_W)),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_B, 
					yvals_lwr=(yvals_lwr_proliferation_B), 
					yvals_upr=(yvals_upr_proliferation_B)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=-wp_mean_B,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_B)),col=specialfill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=(yvals_lwr_clearance_W), 
					yvals_upr=(yvals_upr_clearance_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=10^convert_Ct_logGEML(lod-dp_mean_W),yend=10^convert_Ct_logGEML(lod)),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_B, 
					yvals_lwr=(yvals_lwr_clearance_B), 
					yvals_upr=(yvals_upr_clearance_B)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=0,xend=wr_mean_B,y=10^convert_Ct_logGEML(lod-dp_mean_B),yend=10^convert_Ct_logGEML(lod)),col=specialfill) + 
			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="RNA copies per ml") + 
			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			theme(text=element_text(size=18))
	}

	return(out)

	})

}

# For plotting a normal prior distribution: 
make_normal_prior_df <- function(mean, sd, pmin, pmax, step){
	xmin <- qnorm(pmin, mean=mean, sd=sd)
	xmax <- qnorm(pmax, mean=mean, sd=sd)
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dnorm(xvals, mean=mean, sd=sd)))
	return(out)
}

# For truncating the normal distribution within the bounds specified in the MCMC fit: 
truncnormmean <- function(mu, sigma, T0, T1){
	alpha <- (T0 - mu)/sigma
	beta <- (T1 - mu)/sigma
	out <- mu + 
		(dnorm(alpha, 0, 1) - dnorm(beta, 0, 1))/
		(pnorm(beta, 0, 1) - pnorm(alpha, 0, 1))*sigma
	return(out)
}