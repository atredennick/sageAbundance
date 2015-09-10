##  Script to test STAN poisson model

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
beta <- 0.1
x <- rpois(100, 10)
y <- rpois(100, lambda = exp(alpha+beta*x))
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
  // Likelihood
  y ~ poisson_log(mu);
}
"

####
####  Fit STAN model in chunks of 100 iterations
####
datalist <- list(y=y, x=x, nobs=length(y))
pars <- c("alpha", "beta")

# Compile the model
mcmc_config <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)
mcmc1 <- stan(fit=mcmc_config, data=datalist, pars=pars, chains=1, iter = 2000, warmup = 1000)
traceplot(mcmc1)
print(mcmc1)
