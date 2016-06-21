##  Script for leave-one-year-out cross validation of the
##  sagebrush abundance model. Here we use a simplified version
##  with no spatial random effect to reduce computational overhead.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-09-2015

##  Clear the workspace...
rm(list=ls())


####
####  Set MCMC configurations --------------------------------------------------
####
niters <- 2000
nburns <- 1000
nchains <- 3


####
####  Get leave out year id programmatically -----------------------------------
####
##  For HPC runs
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
year_id <- as.numeric(myargument)

##  For PC tests (comment out for HPC runs!!)
year_id <- 1



####
####  Load libraries -----------------------------------------------------------
####
library(rstan)
library(ggmcmc)



####
####  Read in observation data -------------------------------------------------
####
obs_data <- read.csv("../wy_sagecover_subset_noNA.csv")

# Get data structure right
growD <- subset(obs_data, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Leave out a year
all_years <- unique(growD$Year)
growD <- subset(growD, Year != all_years[year_id])



####
####  Write the STAN model -----------------------------------------------------
####
model_string <- "
data{
  int<lower=0> nobs; // number of observations
  int<lower=0> nyrs; // number of years
  int<lower=0> ncovs; // number of climate covariates
  int<lower=0> yrid[nobs]; // year id
  int y[nobs]; // observation vector
  vector[nobs] lag; // lag cover vector
  matrix[nobs,ncovs] X; // climate covariates
}
parameters{
  real int_mu;
  real<lower=0> beta_mu;
  real<lower=0.000001> sig_yr;
  vector[ncovs] beta;
  vector[nyrs] int_yr;
}
transformed parameters{
  vector[nobs] mu;
  vector[nobs] climEffs;
  climEffs <- X*beta;
  for(n in 1:nobs)
  mu[n] <- int_yr[yrid[n]] + beta_mu*lag[n] + climEffs[n];
}
model{
  // Priors
  sig_yr ~ uniform(0,10);
  int_mu ~ normal(0,100);
  int_yr ~ normal(int_mu, sig_yr);
  beta_mu ~ normal(0,10);
  beta ~ normal(0,10);
  // Likelihood
  y ~ poisson_log(mu);
}
"



####
####  Subset data by removing pixels with 0s -----------------------------------
####
"%w/o%" <- function(x, y) x[!x %in% y] # x without y

# Remove pixels that have 0s in the time series
zids <- growD[which(growD$Cover==0),"ID"]
zids <- unique(c(zids, growD[which(growD$CoverLag==0),"ID"]))
pixels_to_keep <- growD$ID %w/o%  zids
modelD <- growD[which(growD$ID %in% pixels_to_keep),]



####
####  Send data to STAN function for fitting -----------------------------------
####
# Set up new, continuous cell IDs
num_ids <- length(unique(modelD$ID))
num_reps <- nrow(modelD)/num_ids
new_ids <- rep(c(1:num_ids), num_reps)
modelD$newID <- new_ids

y <- modelD$Cover
lag <- log(modelD$CoverLag)
cellid <- modelD$newID
X <- modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X <- scale(X, center = TRUE, scale = TRUE)
yrid <- as.numeric(as.factor(modelD$Year))
nyrs <- length(unique(yrid))

inits <- list()
inits[[1]] <- list(int_mu = 1, beta_mu = 0.5, beta = rep(0, ncol(X)),
                   sig_yr=0.02, int_yr = rep(0,nyrs))
inits[[2]] <- list(int_mu = 2, beta_mu = 1, beta = rep(0.5, ncol(X)),
                   sig_yr=0.1, int_yr = rep(-0.5,nyrs))
inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.8, beta = rep(0.2, ncol(X)),
                   sig_yr=0.15, int_yr = rep(0.5,nyrs))

datalist <- list(y=y, lag=lag, nobs=length(lag), X=X, ncovs=ncol(X), 
                 yrid=yrid, nyrs=nyrs)
pars <- c("int_mu", "beta_mu", "beta", "sig_yr", "int_yr")

# Compile the model
mcmc_config <- stan(model_code=model_string, data=datalist,
                    pars=pars, chains=0)
mcmc <- stan(fit=mcmc_config, data=datalist, pars=pars, chains=nchains, 
              iter = niters, warmup = nburns, init=inits)
longmcmc <- ggs(mcmc)
saveRDS(longmcmc, paste0("mcmc_loyo_year_",all_years[year_id], ".RDS"))
