##  Script to run GLMM fit on subset of sagebrush data

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
rawD <- read.csv(paste(datapath, "WY_SAGECoverData_V2check.csv", sep=""))

# Merge in climate data 
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
# growD <- subset(growD, Year<1994) # subset just 10 years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Load knot data
load("../results/Knot_cell_distances.Rdata")
K <- K.data$K

# Remove NA values
rms <- which(is.na(growD$Cover) == TRUE)
rmsK <- which(is.na(subset(growD, Year==1985)["Cover"])==TRUE)
K <- K[-rmsK,]
growD <- growD[-rms, ]

# Set up new, continuous cell IDs
num_ids <- length(unique(growD$ID))
num_reps <- nrow(growD)/num_ids
new_ids <- rep(c(1:num_ids), num_reps)
growD$newID <- new_ids

####
##  Send data to stan function for fitting
####
y = growD$Cover
lag = growD$CoverLag
cellid = growD$newID
X = growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X = scale(X, center = TRUE, scale = TRUE)

inits <- list()
inits[[1]] <- list(int_mu = 1, beta_mu = 0.05, beta = rep(0, ncol(X)),
                   alpha = rep(0,ncol(K.data$K)), sigma=0.05, sig_a=0.05,
                   sig_mu=0.05, lambda=rep(1, length(y)))
inits[[2]] <- list(int_mu = 2, beta_mu = 0.01, beta = rep(0.5, ncol(X)),
                   alpha = rep(0.5,ncol(K.data$K)), sigma=0.02, sig_a=0.005,
                   sig_mu=0.025, lambda=rep(10, length(y)))
inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.02, beta = rep(0.2, ncol(X)),
                   alpha = rep(0.25,ncol(K.data$K)), sigma=0.04, sig_a=0.025,
                   sig_mu=0.5, lambda=rep(5, length(y)))

model_string <- "
  data{
    int<lower=0> nobs; // number of observations
int<lower=0> ncovs; // number of climate covariates
int<lower=0> nknots; // number of interpolation knots
int<lower=0> ncells; // number of cells
int<lower=0> cellid[nobs]; // cell id
int<lower=0> dK1; // row dim for K
int<lower=0> dK2; // column dim for K
int y[nobs]; // observation vector
int lag[nobs]; // lag cover vector
matrix[dK1,dK2] K; // spatial field matrix
matrix[nobs,ncovs] X; // spatial field matrix
}
parameters{
real int_mu;
real<lower=0> beta_mu;
real<lower=0.000001> sig_a;
real<lower=0.000001> sig_mu;
vector[nknots] alpha;
vector[nobs] lambda;
vector[ncovs] beta;
}
transformed parameters{
vector[ncells] eta;
vector[nobs] mu;
vector[nobs] climEffs;
eta <- K*alpha;
climEffs <- X*beta;
for(n in 1:nobs)
mu[n] <- int_mu + beta_mu*lag[n] + climEffs[n] + eta[cellid[n]];
}
model{
// Priors
alpha ~ normal(0,sig_a);
sig_a ~ cauchy(0, 5);
sig_mu ~ cauchy(0, 5);
int_mu ~ normal(0,100);
beta_mu ~ normal(0,10);
beta ~ normal(0,10);
// Likelihood
lambda ~ normal(mu, sig_mu);
y ~ poisson(exp(lambda));
}
"
datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                   cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K),
                   X=X, ncovs=ncol(X))
pars <- c("int_mu", "beta_mu",  "alpha", "beta", "sig_mu", "sig_a")
  
# Compile the model
mcmc_samples <- stan(model_code=model_string, data=datalist,
                       pars=pars, chains=0)
  
Run parallel chains
rng_seed <- 123
sflist <-
  mclapply(1:3, mc.cores=3,
            function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                            seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                            iter=2000, warmup=1000, init=list(inits[[i]])))
fit <- sflist2stanfit(sflist)
long <- ggs(fit)
saveRDS(mcmc, "stanmcmc_sage.RDS")


