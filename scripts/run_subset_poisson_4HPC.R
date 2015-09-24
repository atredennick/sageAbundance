##  Script to run GLMM fit on subset of sagebrush data
##  Uses STAN to fit the model
##
##  Author:      Andrew Tredennick
##  Email:       atredenn@gmail.com
##  Last Update: 9-2-2015


####
####  Set some file paths, etc.
####
datapath <- ""
knotpath <- ""

####
####  Load required libraries
####
library(rstan)
library(reshape2)


####
####  Set chain ID from command line prompt
####
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
chain_id <- as.numeric(myargument)

####
####  Get data
####
fullD <- read.csv(paste0(datapath,"wy_sagecover_subset_noNA.csv"))
if(length(which(is.na(fullD$Cover))) > 0) stop("data contains NA values")

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Load knot data
load(paste0(knotpath,"Knot_cell_distances_subset.Rdata"))
K <- K.data$K



####
####  Write the STAN model
####
model_string <- "
data{
int<lower=0> nobs; // number of observations
int<lower=0> nyrs; // number of years
int<lower=0> ncovs; // number of climate covariates
int<lower=0> nknots; // number of interpolation knots
int<lower=0> ncells; // number of cells
int<lower=0> cellid[nobs]; // cell id
int<lower=0> yrid[nobs]; // year id
int<lower=0> dK1; // row dim for K
int<lower=0> dK2; // column dim for K
int y[nobs]; // observation vector
vector[nobs] lag; // lag cover vector
matrix[dK1,dK2] K; // spatial field matrix
matrix[nobs,ncovs] X; // spatial field matrix
}
parameters{
real int_mu;
real<lower=0> beta_mu;
real<lower=0.000001> sig_a;
real<lower=0.000001> sig_yr;
vector[nknots] alpha;
vector[ncovs] beta;
vector[nyrs] int_yr;
}
transformed parameters{
vector[ncells] eta;
vector[nobs] mu;
vector[nobs] climEffs;
eta <- K*alpha;
climEffs <- X*beta;
for(n in 1:nobs)
mu[n] <- int_yr[yrid[n]] + beta_mu*lag[n] + climEffs[n] + eta[cellid[n]];
}
model{
// Priors
alpha ~ normal(0,sig_a);
sig_a ~ uniform(0,10);
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
####  Function to extract STAN mcmc
####
get_mcmc <- function(S){
  long <- NULL
  sdf <- as.data.frame(S@sim$samples[[1]])
  sdf$Iteration <- 1:dim(sdf)[1]
  s <- melt(sdf, id.vars = "Iteration")
  colnames(s) <- c("Iteration", "Parameter", "value")
  long <- rbind(long, s)
  return(long)
}



####
####  Subset data by removing pixels with 0s
####
"%w/o%" <- function(x, y) x[!x %in% y] # x without y

# Remove pixels that have 0s in the time series
zids <- growD[which(growD$Cover==0),"ID"]
zids <- unique(c(zids, growD[which(growD$CoverLag==0),"ID"]))
pixels_to_keep <- growD$ID %w/o%  zids
modelD <- growD[which(growD$ID %in% pixels_to_keep),]

# Take out those pixel rows from the K matrix
K <- K.data$K
K <- K[-zids,]

# Check to make sure dimensions match
if(dim(K)[1] != length(unique(modelD$ID))) stop("dimension mismatch between K and data")



####
####  Send data to STAN function for fitting
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
                   alpha = rep(0,ncol(K.data$K)), sig_a=0.05, sig_yr=0.02,
                   int_yr = rep(0,nyrs))
inits[[2]] <- list(int_mu = 2, beta_mu = 1, beta = rep(0.5, ncol(X)),
                   alpha = rep(0.5,ncol(K.data$K)), sig_a=0.005, sig_yr=0.1,
                   int_yr = rep(-0.5,nyrs))
inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.8, beta = rep(0.2, ncol(X)),
                   alpha = rep(0.25,ncol(K.data$K)), sig_a=0.025, sig_yr=0.15,
                   int_yr = rep(0.5,nyrs))

datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                 cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K),
                 X=X, ncovs=ncol(X), yrid=yrid, nyrs=nyrs)
pars <- c("int_mu", "beta_mu",  "alpha", "beta", "sig_a", "sig_yr", "int_yr")
  
# Compile the model
mcmc_config <- stan(model_code=model_string, data=datalist,
                    pars=pars, chains=0)
mcmc1 <- stan(fit=mcmc_config, data=datalist, pars=pars, chains=1, 
              iter = 200, warmup = 100, init=list(inits[[chain_id]]))
long <- get_mcmc(mcmc1)
lastones <- subset(long, Iteration==200)
lastmcmc <- mcmc1
saveRDS(long, paste0("POIS_iterchunk1_chain", chain_id, ".RDS"))

for(i in 2:10){
  betatmp=lastones[grep("beta", lastones$Parameter),]
  beta=as.numeric(unlist(betatmp[which(betatmp$Parameter!="beta_mu"),"value"]))
  newinits <- list(list(int_mu=as.numeric(lastones[which(lastones$Parameter=="int_mu"),"value"]),
                        beta_mu=as.numeric(lastones[which(lastones$Parameter=="beta_mu"),"value"]),
                        beta=beta,
                        alpha=as.numeric(unlist(lastones[grep("alpha", lastones$Parameter),"value"])),
                        sig_a=as.numeric(lastones[which(lastones$Parameter=="sig_a"),"value"]),
                        sig_yr=as.numeric(lastones[which(lastones$Parameter=="sig_yr"),"value"]),
                        int_yr=as.numeric(unlist(lastones[grep("int_yr", lastones$Parameter), "value"]))))
  mcmc <- stan(fit = lastmcmc, data = datalist, pars = pars, chains = 1, 
               iter=200, warmup=100, init = newinits)
  long <- get_mcmc(mcmc)
  lastones <- subset(long, Iteration==200)
  lastmcmc <- mcmc
  saveRDS(long, paste0("POIS_iterchunk",i,"_chain", chain_id, ".RDS"))
}



