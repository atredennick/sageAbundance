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
  sig_a ~ uniform(0,10);
  sig_mu ~ uniform(0,10);
  int_mu ~ normal(0,100);
  beta_mu ~ normal(0,10);
  beta ~ normal(0,10);
  // Likelihood
  lambda ~ normal(mu, sig_mu);
  y ~ poisson(exp(lambda));
}
"

####
##  Send data to stan function for fitting
####
y = growD$Cover
lag = growD$CoverLag
K = K.data$K
cellid = growD$ID
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

datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                 cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K),
                 X=X, ncovs=ncol(X))
pars <- c("int_mu", "beta_mu",  "alpha", "beta")
  
  # Compile the model
  mcmc_samples <- stan(model_code=model_string, data=datalist, init = list(inits[[1]]),
                       pars=pars, chains=1, iter=1000, warmup=500)

outs <- ggs(mcmc_samples)
ggs_traceplot(outs, "alpha")
# saveRDS(outs, file = "small_test_mcmc.RDS")


## Look at spatial field
alpha <- outs[grep("alpha", outs$Parameter),]
alpha_means <- ddply(alpha, .(Parameter), summarise,
                     mean = mean(value))
eta <- K%*%alpha_means$mean
oneyear <- subset(growD, Year==1985)
oneyear$eta <- eta

library(ggplot2)
library(gridExtra)
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=16),text=element_text(size=16),
                legend.text=element_text(size=16))
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(oneyear, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=eta, fill=eta))+
  scale_fill_gradientn(colours=myPalette(100))

ggplot(growD, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=Cover, fill=Cover))+
  facet_wrap("Year")+
  scale_fill_gradientn(colours=myPalette(100))+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))


