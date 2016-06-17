##  Script to simulate a linear model and test regularization
##  via cross-validation using STAN

library(rstan)
library(plyr)
library(ggmcmc)

rm(list=ls())

source("waic_fxns.R")

####  Simulate some data ----
nobs <- 200
alpha <- 10
beta1 <- 2
beta2 <- -5
#Simulate correlated predictor variables
R = matrix(cbind(1,0.80,0.2,0.80,1,0.7,0.2,0.7,1),nrow=3)  #correlation structure for predictors- 3x3, thus 3 beta's
U = t(chol(R))
nvars = dim(U)[1]    #3 beta's
set.seed(1432423)
random.normal=rbind(rep(1,nobs),rnorm(nobs,2,1),rnorm(nobs,-5,1))
X = U %*% random.normal
newX = t(X)
X=newX
X[,2:nvars]=scale(X[,2:nvars],center=TRUE, scale=TRUE)  # scaled to mean 0 and variance 1
# X #Design matrix with intercept and two scaled and centered predictor variables.

y <- rnorm(nobs, X%*%c(alpha,beta1,beta2), 20)


####  Write STAN model ----
model_string <- "
data{
  int<lower=0> nobs; // number observations
  vector[nobs] y; // observation vector
  int<lower=0> ncovs; // number of covariates (+ intercept)
  //vector[ncovs] beta_means;
  real<lower=0> beta_tau;
  matrix[nobs,ncovs] x; // predictor vector
  int<lower=0> npreds; // # of points to predict
  matrix[npreds,ncovs] x_holdout;  // data for which to make predictions
  vector[npreds] y_holdout;  // observations for predictions
}
parameters{
  real b1;
  vector[ncovs] b;
  real<lower=0.000001> tau_model;
}
transformed parameters{
  vector[nobs] mu;
  mu <- b1+x*b;
}
model{
  //Likelihood
  y ~ normal(mu, tau_model);

  //Priors
  b1 ~ normal(0, 1000);
  b ~ normal(0, beta_tau);
}
generated quantities {
  vector[npreds] y_hat; // predictions
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  y_hat <- b1+x_holdout*b;
  for(i in 1:npreds)
    log_lik[i] <- normal_log(y_holdout[i], y_hat[i], tau_model);
}
"

####  Compile STAN model ----
x_train <- X[-c(1:50),2:ncol(X)]
y_train <- y[-c(1:50)]
x_hold <- X[c(1:50),2:ncol(X)]
y_hold <- y[c(1:50)]
nobs_train <- nrow(x_train)
datalist <- list(nobs=nobs_train, y=y_train, x=x_train, ncovs=ncol(x_train),
                 npreds=nrow(x_hold), x_holdout=x_hold, y_holdout=y_hold,
                 beta_means=rep(0,ncol(x_train)), 
                 beta_tau=10)
pars=c("b", "log_lik")
mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=1)

####  Loop over regulators and C-V folds; fit model ----
##  Set up precision range

# lambda.set <- exp(c(-1,0,1,2,5,8))
# lambda.set <- exp(seq(-4,2,length.out=10))
# tau_vec <- 1/lambda.set
# ntaus <- length(tau_vec)
# 
# ##  Set up leave-one-out
# leaveout=sample(rep(1:(nobs/2),2)) #2-fold CV
# nleaveouts <- length(unique(leaveout))


library(parallel) 
library(snowfall) 
library(rlecuyer) 

n.beta=24
s2.beta.vec.2=seq(.25,1.5,,n.beta)^2
K=10
cv.s2.grid=expand.grid(1:n.beta,1:K)
n.grid=dim(cv.s2.grid)[1]
fold.idx.mat=matrix(1:nobs,,K)


iterations <- 200
n.mcmc <- iterations*0.5

# lppd <- rep(NA, nleaveouts)
# lppd2 <- matrix(NA, nrow=nleaveouts, ncol=length(tau_vec))
# betas <- matrix(NA, nrow=nleaveouts, ncol=2)
# betas2 <- array(NA, c(nleaveouts, 2, length(tau_vec)))

cps=detectCores()
sfInit(parallel=TRUE, cpus=cps)
sfExportAll()
sfClusterSetupRNG()

cv.fcn <- function(i){
  library(rstan)
  library()
  k=cv.s2.grid[i,2]
  fold.idx=fold.idx.mat[,k] 
  y_hold=y[fold.idx]
  x_train <- X[-fold.idx,2:ncol(X)]
  y_train <- y[-fold.idx]
  x_hold <- X[fold.idx,2:ncol(X)]
  y_hold <- y[fold.idx]
  nobs_train <- nrow(x_train)
  datalist <- list(nobs=nobs_train, y=y_train, x=x_train, ncovs=ncol(x_train),
                   npreds=nrow(x_hold), x_holdout=x_hold, y_holdout=y_hold,
                   beta_means=rep(0,ncol(x_train)), 
                   beta_tau=s2.beta.vec.2[cv.s2.grid[i,1]])
  pars <- c("log_lik")
  fit <- stan(fit=mcmc_samples, data=datalist,
              pars=pars, chains=1, iter = iterations, warmup = n.mcmc)
#   long <- ggs(fit)
#   postmeans <- ddply(long, .(Parameter), summarise,
#                      postmean = mean(value))
#   tmpbetas <- postmeans[grep("b", postmeans$Parameter), "postmean"]
#   betas[j,] <- tmpbetas
  waic_metrics <- waic(fit)
  lppd <- waic_metrics[["elpd_loo"]]
  lppd
}

sfExport("y","X","n.mcmc","s2.beta.vec.2","cv.s2.grid","cv.fcn","fold.idx.mat","mcmc_samples")
tmp.time=Sys.time()
score.list=sfClusterApplySR(1:n.grid,cv.fcn,perUpdate=1)
time.2=Sys.time()-tmp.time
sfStop()

time.2
score.cv.mat=matrix(unlist(score.list),n.beta,K)


# 
# for(i in 1:ntaus){
#   for(j in 1:nleaveouts){
#     temp <- which(leaveout==j)
#     x_train <- X[-temp,2:ncol(X)]
#     y_train <- y[-temp]
#     x_hold <- X[temp,2:ncol(X)]
#     y_hold <- y[temp]
#     nobs_train <- nrow(x_train)
#     
#     datalist <- list(nobs=nobs_train, y=y_train, x=x_train, ncovs=ncol(x_train),
#                      npreds=nrow(x_hold), x_holdout=x_hold, y_holdout=y_hold,
#                      beta_means=rep(0,ncol(x_train)), 
#                      beta_tau=rep(tau_vec[i],ncol(x_train)))
#     pars <- c("b", "log_lik")
#     fit <- stan(fit=mcmc_samples, data=datalist,
#                 pars=pars, chains=1, iter = iterations, warmup = n.mcmc)
#     long <- ggs(fit)
#     postmeans <- ddply(long, .(Parameter), summarise,
#                        postmean = mean(value))
#     tmpbetas <- postmeans[grep("b", postmeans$Parameter), "postmean"]
#     betas[j,] <- tmpbetas
#     waic_metrics <- waic(fit)
#     lppd[j] <- waic_metrics[["elpd_loo"]]
#     # lppd[j] <- sum(postmeans[grep("lppd", postmeans$Parameter), "postmean"])
#   }
#   lppd2[,i] <- lppd
#   betas2[,,i] <- betas
# }

## Sum across C-Vs
temp1 <- colSums(lppd2)
opt.reg <- which(temp1==max(temp1))
plot(log(1/tau_vec), temp1, type="l")
abline(v=log(1/tau_vec[opt.reg]), col="grey")
temp2 <- colMeans(betas2)
matplot(t(temp2), type="l")
abline(h=0, lty=2)
abline(v=opt.reg, col="grey")


### Fit optimal predictive model

# model_string <- "
# data{
#   int<lower=0> nobs; // number observations
#   vector[nobs] y; // observation vector
#   int<lower=0> ncovs; // number of covariates (+ intercept)
#   vector[ncovs] beta_means;
#   vector[ncovs] beta_tau;
#   matrix[nobs,ncovs] x; // predictor vector
# }
# parameters{
#   real b1;
#   vector[ncovs] alpha;
#   //vector[ncovs] b;
#   real<lower=0> tau_model;
#   corr_matrix[ncovs] Tau;
# }
# transformed parameters{
#   matrix[ncovs, ncovs] L;
#   vector[nobs] mu;
#   vector[ncovs] b;
#   L <- cholesky_decompose(Tau);
#   b <- beta_means + beta_tau .* (L * alpha);
#   //matrix[ncovs,ncovs] Sigma;
#   mu <- b1+x*b;
#   //Sigma <- diag_post_multiply(Omega, beta_tau);
# }
# model{
#   //Likelihood
#   y ~ normal(mu, tau_model);
# 
#   //Priors
#   //Omega ~ lkj_corr(1);
#   alpha~normal(0,1);
#   b1 ~ normal(0, 1000);
#   //b ~ multi_normal(beta_means, Sigma);
# }
# "
# # tau <- tau_vec[opt.reg]
# tau <- 1
# x_train <- X[,2:ncol(X)]
# y_train <- y
# nobs_train <- nrow(x_train)
# datalist <- list(nobs=nobs_train, y=y_train, x=x_train, ncovs=ncol(x_train),
#                  beta_means=rep(0,ncol(x_train)), 
#                  beta_tau=rep(tau,ncol(x_train)))
# pars=c("b", "b1", "Tau", "L", "alpha")
# mcmc_samples <- stan(model_code=model_string, data=datalist,
#                      pars=pars, chains=1)
# plot(mcmc_samples)
# print(mcmc_samples)
