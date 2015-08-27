##  Script to test dividing up STAN MCMC into chunks

# Clear workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(rstan)
library(ggmcmc)

####
####  Simulate some data
####
alpha <- 10
beta <- 5
sigma <- 2
x <- rnorm(100,0,1)
y <- rnorm(100,alpha+beta*x, sigma)
plot(x, y)

####
####  Write STAN model
####
model_string <- "
data{
  int<lower=0> nobs; // number of observations
  vector[nobs] y; // observation vector
  vector[nobs] x; // covariate vector
}
parameters{
  real alpha;
  real beta;
  real<lower=0.000001> sigma;
}
transformed parameters{
  vector[nobs] mu;
  mu <- alpha + beta*x;
}
model{
  // Priors
  alpha ~ normal(0,1000);
  beta ~ normal(0,1000);
  sigma ~ uniform(0,100);
  // Likelihood
  y ~ normal(mu, sigma);
}
"

####
####  Fit STAN model in chunks of 100 iterations
####
datalist <- list(y=y, x=x, nobs=length(y))
pars <- c("alpha", "beta", "sigma")

# Compile the model
mcmc_config <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)
mcmc1 <- stan(fit=mcmc_config, data=datalist, pars=pars, chains=1, iter = 200, warmup = 100)
traceplot(mcmc1)
long <- ggs(mcmc1)
lastones <- subset(long, Iteration==100)

newinits <- list(list(alpha=as.numeric(lastones[which(lastones$Parameter=="alpha"),"value"]),
                 beta=as.numeric(lastones[which(lastones$Parameter=="beta"),"value"]),
                 sigma=as.numeric(lastones[which(lastones$Parameter=="sigma"),"value"])))
mcmc2 <- stan(fit = mcmc1, data = datalist, pars = pars, chains = 1, iter=200, warmup=100, init = newinits)
traceplot(mcmc2)
long2 <- ggs(mcmc2, inc_warmup = TRUE)
firstones <- subset(long2, Iteration==1)

firstones
lastones

####
####  Combine MCMC chains and plot
####
long2$Iteration <- long2$Iteration+100
long_mcmc <- rbind(long, long2)
ggs_traceplot(long_mcmc)

