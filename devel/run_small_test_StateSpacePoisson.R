##  Script to run GLMM fit on subset of sagebrush data

# Clear workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(sageAbundance)
library(rstan)
library(ggmcmc)
library(parallel)
library(reshape2)
library(plyr)
library(gridExtra)

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
####  Write the STAN model
####
model_string <- "
data{
  int<lower=0> nyrs; // number of observation years
  int<lower=0> nknots; // number of interpolation knots
  int<lower=0> ncells; // number of cells
  int<lower=0> dK1; // row dim for K
  int<lower=0> dK2; // column dim for K
  int y[nyrs,ncells]; // observation matrix
  matrix[dK1,dK2] K; // spatial field matrix
  vector[nyrs] X; // climate covariate matrix
  vector[ncells] startvec;
}
parameters{
  real int_mu;
  real<lower=0> beta_mu;
  real<lower=0.0001> sig_a;
  vector[nknots] alpha;
  real beta1;
}
transformed parameters{
  vector[ncells] eta;
  matrix[nyrs, ncells] mu;
  eta <- K*alpha;
  for(i in 1:ncells)
    mu[1,i] <- startvec[i];
  for(t in 2:nyrs){
    for(i in 1:ncells){
      mu[t,i] <- int_mu + beta_mu*mu[t-1,i] + X[t]*beta1 + eta[i];
    }
  }
}
model{
  // Priors
  alpha ~ normal(0,sig_a);
  int_mu ~ normal(0,100);
  beta_mu ~ normal(0,10);
  beta1 ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  // Likelihood
  for(t in 1:nyrs)
    for(i in 1:ncells){
      y[t,i] ~ poisson_log(mu[t,i]);
    }
    
}
"

####
####  Send data to STAN function for fitting
####
wide_data <- dcast(growD, Year~ID, value.var = "Cover")
y <- wide_data[,2:ncol(wide_data)]
startvec <- as.numeric(log(y[1,]))
K <- K.data$K
cellid <- growD$ID
clim_covs <- c("pptLag")
X <- climD[which(climD$year %in% unique(growD$Year)),clim_covs]
X <- as.numeric(scale(X, center = TRUE, scale = TRUE))


datalist <- list(y=y, nyrs=nrow(y), ncells=length(unique(cellid)),
                 cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K),
                 X=X, startvec=startvec)
pars <- c("int_mu",  "alpha", "beta1", "beta_mu")
  
# Compile the model
mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=0)
test <- stan(fit=mcmc_samples, data=datalist, pars=pars,
             chains=1, iter=2000, warmup=1000)
print(test)
traceplot(test)

fit <- ggs(test)
