##  Script to run sage abundance model through time, starting
##  at the last observation. Requires time series of climate
##  projections.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-08-2015


### Clear the workspace
rm(list=ls())



####
####  Set some global simulation settings --------------------------------------
####
parameter_reps <- 100
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")


####
####  Load necessary libraries -------------------------------------------------
####
library(sageAbundance)
library(plyr)
library(reshape2)

#TODO: get sage Abundance functions to load!! Until then...
source("../R/clim_proj_format_funcs.R")



####
####  Read in percent cover file, climate projections, MCMC output, K.data -----
####
obs_data <- read.csv("../data/wy_sagecover_subset_noNA.csv")
obs_data <- subset(obs_data, Year>1984) # subsets out NA CoverLag values
temp_projs <- readRDS("../data/CMIP5_yearly_project_temperature.RDS")
ppt_projs <- readRDS("../data/CMIP5_yearly_project_precipitation.RDS")
mcmc_outs <- readRDS("../results/poissonSage_randYear_mcmc.RDS")
load("../results/Knot_cell_distances_subset.Rdata")
K <- K.data$K



####
####  Subset out observed climate; get scaling mean and sd ------------------------
####
obs_clim <- obs_data[,c("Year",clim_vars)]
obs_clim <- unique(obs_clim)
obs_clim_means <- colMeans(obs_clim[,clim_vars])
obs_clim_sds <- apply(obs_clim[,clim_vars], 2, sd)
obs_clim_scalers <- data.frame(variable = clim_vars,
                               means = obs_clim_means,
                               sdevs = obs_clim_sds)

####
####  Define simulation function -----------------------------------------------
####
iterate_sage <- function(N, int, beta.dens, beta.clim, eta, weather){
  dens.dep <- beta.dens*log(N)
  clim.effs <- sum(beta.clim*weather)
  Nout <- exp(int + dens.dep + clim.effs + eta)
  return(Nout)
}



####
####  Begin simulation set up --------------------------------------------------
####
last_year <- max(obs_data$Year)
last_obs <- subset(obs_data, Year == last_year)
all_models <- unique(temp_projs$model)
all_scenarios <- unique(temp_projs$scenario)
sim_years <- c((last_year+1):max(temp_projs$Year))
num_sims <- length(sim_years)

nchains <- length(unique(mcmc_outs$chain))
niters <- length(unique(mcmc_outs$mcmc_iter))

####
####  Begin looping: parameters within years within scenario within model ------
####
for(do_model in all_models){
  
  for(do_scenario in all_scenarios){
    
    temp_now <- subset(temp_projs, scenario==do_scenario & model==do_model)
    ppt_now <- subset(ppt_projs, scenario==do_scenario & model==do_model)
    climate_now <- format_climate(tdata = temp_now, 
                                  pdata = ppt_now, 
                                  years = sim_years)
    
    # Scale climate predictors
    climate_now["pptLag"] <- (climate_now["pptLag"] - obs_clim_means["pptLag"])/obs_clim_sds["pptLag"]
    climate_now["ppt1"] <- (climate_now["ppt1"] - obs_clim_means["ppt1"])/obs_clim_sds["ppt1"]
    climate_now["ppt2"] <- (climate_now["ppt2"] - obs_clim_means["ppt2"])/obs_clim_sds["ppt2"]
    climate_now["TmeanSpr1"] <- (climate_now["TmeanSpr1"] - obs_clim_means["TmeanSpr1"])/obs_clim_sds["TmeanSpr1"]
    climate_now["TmeanSpr2"] <- (climate_now["TmeanSpr2"] - obs_clim_means["TmeanSpr2"])/obs_clim_sds["TmeanSpr2"]
    
    # Create storage matrix for population
    n_save <- array(dim = c(num_sims+1, parameter_reps, nrow(last_obs)))
    n_save[1,,] <- last_obs$Cover # set first record to last observation
    
    for(t in 2:num_sims){
      
      weather <- climate_now[t-1, clim_vars] #t-1 since clim data starts in 2012
      
      for(i in 1:parameter_reps){
        randchain <- sample(1:nchains, 1)
        randiter <- sample(1:niters, 1)
        params_now <- subset(mcmc_outs, chain==randchain & mcmc_iter==randiter)
        alphas <- params_now[grep("alpha", params_now$Parameter), "value"]
        int_mu <- params_now[grep("int_mu", params_now$Parameter), "value"]
        beta_mu <- params_now[grep("beta_mu", params_now$Parameter), "value"]
        beta_clims <- params_now[grep("beta.", params_now$Parameter), "value"][2:6]
        eta_now <-  K%*%as.numeric(unlist(alphas))
        
        n_save[t,i,] <- iterate_sage(N = n_save[t-1,i,], 
                                     int = int_mu, 
                                     beta.dens = beta_mu,
                                     beta.clim = beta_clims,
                                     eta = eta_now,
                                     weather = weather)
        
      } # next parameter set
      
    } # next year (t)
    
    file <- paste0(do_model,"_", do_scenario, "_yearly_forecasts.RDS")
    path <- "../results/yearlyforecasts/"
    saveRDS(n_save, paste0(path,file))
    
  } # next scenario
  
} # next model




