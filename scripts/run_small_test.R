##  Script to run GLMM fit on subset of sagebrush data

# Clear workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(sageAbundance)
library(rstan)
library(ggmcmc)
library(parallel)

####
##  Get data
####
datapath <-  "/Users/atredenn/Dropbox/sageAbundance_data/"
climD <- read.csv(paste(datapath,
                        "/studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv",
                        sep=""))
rawD <- read.csv(paste(datapath, "WY_SAGECoverData_V2small.csv", sep=""))

# Merge in climate data 
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Load knot data
load("../results/Knot_cell_distances_smallSet.Rdata")


####
## Quick glm for initial values
####
# mod <- glm(Cover ~ CoverLag, family="poisson", data=growD)
# summary(mod)

model_string <- "
data{
  int<lower=0> nobs; // number of observations
  int<lower=0> nknots; // number of interpolation knots
  int<lower=0> ncells; // number of cells
  int<lower=0> cellid[nobs]; // cell id
  int<lower=0> dK1; // row dim for K
  int<lower=0> dK2; // column dim for K
  int y[nobs]; // observation vector
  int lag[nobs]; // lag cover vector
  matrix[dK1,dK2] K; // spatial field matrix
}
parameters{
  real int_mu;
  real<lower=0> beta_mu;
  real<lower=0> sig_a;
  real<lower=0> sig_mu;
  vector[nknots] alpha;
  vector[nobs] lambda;
}
transformed parameters{
  vector[ncells] eta;
  vector[nobs] mu;
  eta <- K*alpha;
  for(n in 1:nobs)
    mu[n] <- int_mu + beta_mu*lag[n] + eta[cellid[n]];
}
model{
  // Priors
  alpha ~ normal(0,sig_a);
  sig_a ~ uniform(0,10);
  sig_mu ~ uniform(0,10);
  int_mu ~ normal(0,100);
  beta_mu ~ normal(0,10);
  // Likelihood
  lambda ~ normal(mu, sig_mu);
  y ~ poisson(exp(lambda));
}
"

####
##  Send data to stan function for fitting
####
inits <- list()
inits[[1]] <- list(int_mu = 1, beta_mu = 0.05, 
                   alpha = rep(0,ncol(K.data$K)), sigma=0.05, sig_a=0.05)
inits[[2]] <- list(int_mu = 2, beta_mu = 0.01, 
                   alpha = rep(0.5,ncol(K.data$K)), sigma=0.02, sig_a=0.005)
inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.02, 
                   alpha = rep(0.25,ncol(K.data$K)), sigma=0.04, sig_a=0.025)

y = growD$Cover
lag = growD$CoverLag
K = K.data$K
cellid = growD$ID

datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                   cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K))
pars <- c("int_mu", "beta_mu",  "alpha")
  
  # Compile the model
  mcmc_samples <- stan(model_code=model_string, data=datalist,
                       pars=pars, chains=3, iter=1000, warmup=500)

outs <- ggs(mcmc_samples)
saveRDS(outs, file = "small_test_mcmc.RDS")










# 
# 
# 
# inits[[2]] <- list(int_mu = 2, beta_mu = 0.005, 
#                    alpha = rep(0.003,ncol(K.data$K)), sigma=0.05)
# mcmc <- model_nocovars(y = growD$Cover, lag = growD$CoverLag, K = K.data$K, 
#                        cellid = growD$ID, iters = 100, warmup = 25, 
#                        nchains = 2, inits=inits)
# ggs_traceplot(mcmc, "int")
# 
# 
