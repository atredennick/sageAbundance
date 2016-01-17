##  Sourcing script for analysis and results 
##  from Tredennick et al. 201x
##
##  Date created: 1-17-2016
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com

##  NOTICE: if you wish to reproduce our results, take care to change directory
##  paths to match your local machine. The easiest way to do this is to clone
##  the repository from GitHub (link below) and just change directory paths
##  in this sourcing script. We provide the MCMC output from the fitting stage,
##  so no need to re-fit the model unless you just want to. Same goes for
##  subsetting the data and calculating the expansion matrix K. So, we comment
##  those out here.


### Subset data from full 5x5km raster to 5x2.5km working extent and
### place and parameterize knots (basis function matrix, K)
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
# source("subsetData_getKnots.R")


### Fit sage abundance GLMM
### This will take a very, VERY long time on a PC; 
### we used Utah State University's High Performance Computing System. 
### Users should NOT attempt to fit this model on a PC.
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
# source("run_subset_poisson_4HPC.R")


### Make study area figure (Figure 1)
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
source("make_study_map.R")


### Calculate climate changes and plot climate time series 
### (Table 1 and Figure 2)
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
source("plot_projected_weather.R")


### Plot posterior climate covariates (Figure 3);
### Run equilibrium simulations and plot (Figure 4);
### Calculate in-sample RMSE and correlation;
### Make Appendix plots of eta 
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
source("ms_scripts.R")


### Run and plot climate change equilibrium simulations (Figure 5)
### This will take awhile...hour or more on PC
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
source("sim_sage_byYear_byModel.R") # simulations and Figure 5


### Run temporally-explicit forecasts (Figure 6)
### This will take awhile...hour or more on PC
setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
source("sim_sage_byYear_byModel.R") # simulations
source("plot_temporal_forecasts.R") # plot Figure 6


