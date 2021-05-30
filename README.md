# Data and code associated with "Densely sampled viral trajectories for SARS-CoV-2 variants B.1.1.7 and B.1.429"
S.M. Kissler`*`, J.R. Fauver`*`, C. Mack`*`, C. Tai, M. Breban, A. Watkins, R.M. Samant, D.J. Anderson, D.D. Ho, N. Grubaugh`+`, Y.H. Grad`+`

`*` denotes equal contribution

`+` denotes co-senior author

Correspondence: skissler@hsph.harvard.edu

__run_analysis.R__ is the main analysis file. It calls all other functions in the proper order. It calls data from the data/ directory, saves key figures into the figures/ directory, and saves the fitted distributions into the output/ directory. 

__utils.R__ contains key functions needed for the analysis. 

__utils_private.R__ contains definitions and functions to be hidden from Github. Currently, this just contains one line specifying the path to an external directory (stored as variable 'extdrive') where the MCMC output is stored.  

__set_global_pars.R__ defines parameters needed throughout the code (currently just the PCR limit of detection). 

__ct_dat_refined.csv__ is the main data file. It contains a person ID, the status (Yes/No) for all current variants of interest/variants of concern, the Test Date Index (days from that person's minimum recorded Ct value), the Roche cobas Target 1 Ct value adjusted to the Yale lab scale, an indicator for whether or not the individual reported symptoms (1/0 - note that the timing of symptoms was not recorded), the lineage, the individual's age group, and a unique observation ID.

__set_run_pars.R__ specifies a list of input parameters for running different versions of the model (e.g., which rows we're contrasting agains the others, which rows to exclude, which priors to use).

__fit_posteriors_preamble.R__ defines key data structures for fitting the viral trajectories.

__fit_posteriors.R__ runs the MCMC sampling scheme to fit the viral trajectory parameters.

__fit_posteriors_special.stan__ contains the Stan program for fitting the viral trajectory parameters. 

__extract_parameters.R__ extracts MCMC draws from the Stan fit into a data table. 

__make_figures.R__ generates the figures from the manuscript. 

__summarise_dists.R__ summarizes the posterior distributions from the MCMC fit in a more convenient data structure. 

__save_figures.R__ writes the figures to .pdf or .png. The directory for saving the figures is speified in __run_analysis.R__ (default is the figures/ directory). 

__diagnose_fit.R__ stores key diagnostics from the MCMC fit, including the Gelman Rhat statistic and the number of divergent transitions. 

__make_figures_overall.R__ generates figures comparing fits from across runs (currently set to compare B.1.1.7 trajectories and B.1.429 trajectories against non-VOI/VOC trajectories). 


All code is available under the GNU General Public License, version 3, included in this repository under the following terms: 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

