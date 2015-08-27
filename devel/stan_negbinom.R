##  Script to test STAN negative binomial model

# Clear workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(rstan)
library(ggmcmc)

####
####  Look at difference between poisson and neg binomial
####
plot(density(rpois(10000, 10), adjust=2), col="red")
lines(density(rnbinom(10000, size = 5, mu = 10), adjust=2))

####
####  Simulate some data
####
alpha <- 1
beta <- 0.07
phi <- 1.5
x <- rpois(100,10)
y <- rnbinom(100, mu = exp(alpha+beta*x), size = phi)
plot(x, y)
abline(0,1)

####
####  Write STAN model
####
model_string <- "
data{
  int<lower=0> nobs; // number of observations
  int y[nobs]; // observation vector
  int x[nobs]; // covariate vector
}
parameters{
  real alpha;
  real beta;
  real<lower=0.000001> phi;
}
transformed parameters{
  vector[nobs] mu;
  for(n in 1:nobs)
    mu[n] <- alpha + beta*x[n];
}
model{
  // Priors
  alpha ~ normal(0,1000);
  beta ~ normal(0,1000);
  phi ~ cauchy(0,3);
  // Likelihood
  y ~ neg_binomial_2_log(mu, phi);
}
"

####
####  Fit STAN model in chunks of 100 iterations
####
datalist <- list(y=y, x=x, nobs=length(y))
pars <- c("alpha", "beta", "phi")

# Compile the model
mcmc_config <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)
mcmc1 <- stan(fit=mcmc_config, data=datalist, pars=pars, chains=3, iter = 2000, warmup = 1000)
traceplot(mcmc1)
print(mcmc1)
