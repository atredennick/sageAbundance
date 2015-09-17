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
  int<lower=0> nobs; // number of observations
  int<lower=0> npres; // number of observations
  int<lower=0> nabs; // number of observations
  int<lower=0> pres;
  int<lower=0> absd;
  int<lower=0> ncovs; // number of climate covariates
  int<lower=0> nknots; // number of interpolation knots
  int<lower=0> ncells; // number of cells
  int<lower=0> cellid[nobs]; // cell id
  int<lower=0> dK1; // row dim for K
  int<lower=0> dK2; // column dim for K
  real y[nobs]; // observation vector
  //real lag[nobs]; // lag cover vector
  matrix[dK1,dK2] K; // spatial field matrix
  matrix[nobs,ncovs] X; // spatial field matrix
}
parameters{
  real int_mu;
  real a;
  real b;
  //real<lower=0> beta_mu;
  real<lower=0.000001> sig_a;
  real<lower=0.000001> sig_mu;
  vector[nknots] alpha;
  vector[ncovs] beta;
}
transformed parameters{
  vector[ncells] eta;
  vector[nobs] mu;
  vector[nobs] climEffs;
  vector[npres] Y_a;
  vector<lower=0, upper=1>[nabs] proba;            // parameter for beta distn
  int<lower=0> npatches[npres];         // parameter for beta distn
  eta <- K*alpha;
  climEffs <- X*beta;
  for(n in 1:nobs)
    mu[n] <- exp(int_mu + climEffs[n] + eta[cellid[n]]);
}
model{
  // Priors
  alpha ~ normal(0,sig_a);
  sig_a ~ uniform(0,10);
  sig_mu ~ uniform(0,10);
  int_mu ~ normal(0,100);
  //beta_mu ~ normal(0,10);
  beta ~ normal(0,10);
  a ~ gamma(2,5);
  b ~ gamma(2,5);

  // Likelihoods
    //Strictly positive percent cover values
    for(k in 1:npres){
      npatches[k] ~ poisson(mu[pres[k]]);
      Y_a[k] <- -a*npatches[k];
      Y[pres[k]] ~ gamma(Y_a[k], b);
    }

    //Evaluation of probabilities of zeros
    for(j in 1:nabs){
      proba[j] <- 1 - exp(-mu[absd[j]]);
      Y[abs[j]] ~ bernoulli(proba[j]);
    }
    
}
"

####
####  Send data to STAN function for fitting
####
y = growD$Cover
lag = growD$CoverLag
K = K.data$K
cellid = growD$ID
X = growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X = scale(X, center = TRUE, scale = TRUE)
abs <- which(y==0)
pres <- which(y>0)
npres <- length(pres)
nabs <- length(abs)


# inits <- list()
# inits[[1]] <- list(int_mu = 1, beta_mu = 0.05, beta = rep(0, ncol(X)),
#                    alpha = rep(0,ncol(K.data$K)), sigma=0.05, sig_a=0.05,
#                    sig_mu=0.05, lambda=rep(1, length(y)), phi=10)
# inits[[2]] <- list(int_mu = 2, beta_mu = 0.01, beta = rep(0.5, ncol(X)),
#                    alpha = rep(0.5,ncol(K.data$K)), sigma=0.02, sig_a=0.005,
#                    sig_mu=0.025, lambda=rep(10, length(y)), phi=20)
# inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.02, beta = rep(0.2, ncol(X)),
#                    alpha = rep(0.25,ncol(K.data$K)), sigma=0.04, sig_a=0.025,
#                    sig_mu=0.5, lambda=rep(5, length(y)), phi=100)

datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                 cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K),
                 X=X, ncovs=ncol(X), nabs=nabs, npres=npres, pres=pres, absd=abs)
pars <- c("int_mu",  "alpha", "beta", "phi")
  
# Compile the model
mcmc_samples <- stan(model_code=model_string, data=datalist, init = list(inits[[1]]),
                     pars=pars, chains=0)
# test <- stan(fit=mcmc_samples, data=datalist, pars=pars,
#              chains=1, iter=200, warmup=100)

