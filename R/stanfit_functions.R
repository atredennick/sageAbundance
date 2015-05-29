#' Estimate model parameters in GLMM using rstan with no environmental covariates
#' 
#' @author Andrew Tredennick
#' @param y Vector of responses (cover observations)
#' @param lag Vector of lag response for temporal response (previous year's cover)
#' @param K Spatial field matrix
#' @param cellid Vector of unique ID's for each spatial cell, replicate through time
#' @param iters Number of MCMC iterations to run
#' @param inits A list of lists whose length is equal to number of chains. Elements
#'              include initial values for parameters to be estimated.
#' @param warmup Number of MCMC iterations to throw out before sampling (< iters)
#' @param nchains Number of MCMC chains to sample
#' @export
#' @return A list with ggs object with MCMC values from rstan and 
#'          Rhat values for each parameter.
model_nocovars <- function(y, lag, K, cellid, iters=2000, warmup=1000, nchains=1, inits){
  ##  Stan C++ model
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
    real beta_mu;
    vector[nknots] alpha;
  }
  transformed parameters{
    vector[ncells] eta;
    vector[nobs] lambda;
    eta <- K*alpha;
    for(n in 1:nobs)
      lambda[n] <- exp(int_mu + beta_mu*lag[n] + eta[cellid[n]]);
  }
  model{
    // Priors
    for(k in 1:nknots){
      alpha[k] ~ normal(0,1000);
    }
    int_mu ~ normal(0,1000);
    beta_mu ~ normal(0,1000);
    // Likelihood
    y ~ poisson(lambda);
  }
  "
  datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                   cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K))
  pars <- c("int_mu", "beta_mu",  "alpha")
  fit <- stan(model_code=model_string, data=datalist, iter=iters,
              warmup=warmup, pars=pars, chains=nchains, init=inits)
  long <- ggs(fit)
#   stansumm <- as.data.frame(summary(fit)["summary"])
#   rhats <- stansumm["summary.Rhat"]
#   return_list <- list(mcmc=long, rhat=rhats)
  return(long)
}