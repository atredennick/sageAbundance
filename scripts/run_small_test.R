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
  vector[ncovs] beta;
  real<lower=0.000001> phi; // neg. binomial dispersion parameter
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
  phi ~ cauchy(0, 3);
  // Likelihood
  y ~ neg_binomial_2_log(mu, phi);
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

inits <- list()
inits[[1]] <- list(int_mu = 1, beta_mu = 0.05, beta = rep(0, ncol(X)),
                   alpha = rep(0,ncol(K.data$K)), sigma=0.05, sig_a=0.05,
                   sig_mu=0.05, lambda=rep(1, length(y)), phi=10)
inits[[2]] <- list(int_mu = 2, beta_mu = 0.01, beta = rep(0.5, ncol(X)),
                   alpha = rep(0.5,ncol(K.data$K)), sigma=0.02, sig_a=0.005,
                   sig_mu=0.025, lambda=rep(10, length(y)), phi=20)
inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.02, beta = rep(0.2, ncol(X)),
                   alpha = rep(0.25,ncol(K.data$K)), sigma=0.04, sig_a=0.025,
                   sig_mu=0.5, lambda=rep(5, length(y)), phi=100)

datalist <- list(y=y, lag=lag, nobs=length(lag), ncells=length(unique(cellid)),
                 cellid=cellid, nknots=ncol(K), K=K, dK1=nrow(K), dK2=ncol(K),
                 X=X, ncovs=ncol(X))
pars <- c("int_mu", "beta_mu",  "alpha", "beta", "phi")
  
# Compile the model
mcmc_samples <- stan(model_code=model_string, data=datalist, init = list(inits[[1]]),
                     pars=pars, chains=0)

# Fit model in parallel
rng_seed <- 12345
sflist <-
  mclapply(1:3, mc.cores=3,
           function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                            seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                            iter=2000, warmup=1000, init=list(inits[[i]]), thin=5))
fit <- sflist2stanfit(sflist)
outs <- ggs(fit)
# ggs_traceplot(outs, "alpha")
# ggs_traceplot(outs, "beta")
saveRDS(outs, file = "../results/small_test_mcmc.RDS")


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
# ggplot(oneyear, aes(x=Lon, y=Lat))+
#   geom_raster(aes(z=eta, fill=eta))+
#   scale_fill_gradientn(colours=myPalette(100))
# 
# ggplot(growD, aes(x=Lon, y=Lat))+
#   geom_raster(aes(z=Cover, fill=Cover))+
#   facet_wrap("Year")+
#   scale_fill_gradientn(colours=myPalette(100))+
#   tmp.theme+
#   theme(strip.background=element_rect(fill="white"))


####
####  Run population simulation
####
mean_params <- ddply(outs, .(Parameter), summarise,
                     mean(value))
alphas <- mean_params[grep("alpha", mean_params$Parameter),"..1"]
betas <- mean_params[grep("beta", mean_params$Parameter),"..1"][2:6]
eta <- K%*%alphas
time.steps <- 1000
pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
clim_sim <- climD[climD$year %in% unique(growD$Year),]
X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- exp(mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,]) + 
    sum(betas*Xtmp)
  tmp.mu <- tmp.mu + eta
  ex.mat[t,] <- rnbinom(ncol(ex.mat), mu=tmp.mu, size = mean_params[mean_params$Parameter=="phi","..1"])
}

# matplot(c(1:time.steps), ex.mat, type="l", col="grey")
# hist(y, freq = FALSE)
# lines(density(ex.mat, adjust = 4), col="red", lwd=2)
# mean(y)
# mean(ex.mat)


####
####  Plot spatial equilibrium
####
obs.equil <- ddply(growD, .(Lon, Lat), summarise,
                   cover = mean(Cover))
obs.equil <- obs.equil[with(obs.equil, order(-Lat, Lon)), ]
obs.equil[which(obs.equil[,"cover"]>10), "cover"] <- NA
eq.pixels <- colMeans(ex.mat)
sp.equil <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
                       Lat=subset(growD, Year==1985)$Lat,
                       pred.cover=eq.pixels,
                       obs.cover=obs.equil$cover)
colnames(sp.equil)[3:4] <- c("Predicted", "Observed") 
sp.equil2 <- melt(sp.equil, id.vars = c("Lon", "Lat"))



png("../results/obs_pred_spatial_small.png", width = 5, height=5, units = "in", res=200)
g1 <- ggplot(sp.equil2, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=value, fill=value))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover")+
  facet_wrap("variable")+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))

g2 <- ggplot(sp.equil)+
  geom_histogram(aes(x=Observed, y = ..density..), col="white", fill="grey45", binwidth=1)+
  geom_line(stat="density", aes(x=Predicted, y= ..density..), 
            col="black", size=1.5, adjust=15, linetype=2)+
  xlab("Percent Cover")+
  ylab("Estimated Kernel Density")+
  theme_bw()
gout <- grid.arrange(g1, g2, ncol=1, nrow=2)
dev.off()

####
####  Plot estimated coefficient densities
####
betas <- outs[grep("beta", outs$Parameter),]
betas <- as.data.frame(subset(betas, Parameter!="beta_mu"))
beta.names <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")

png("../results/climeff_dens_small.png", width = 6, height=3.5, units = "in", res=200)
ggplot(betas)+
  geom_line(stat="density", aes(x=value, y=..density..,color=Parameter), size=1.5)+
  scale_color_discrete(labels=beta.names, name="Coefficient")+
  geom_vline(xintercept=0, linetype=2)+
  xlab("Standardized Coefficient Value")+
  ylab("Posterior Density")+
  theme_bw()
dev.off()




####
####  Run climate simulations
####
projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
                  "deltaPpt"=c(1.0894,1.0864,1.110),
                  "deltaTspr"=c(2.975,3.134,4.786))
mean_params <- ddply(outs, .(Parameter), summarise,
                     mean(value))
alphas <- mean_params[grep("alpha", mean_params$Parameter),"..1"]
betas <- mean_params[grep("beta", mean_params$Parameter),"..1"][2:6]
eta <- K%*%alphas
time.steps <- 1000

pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[1,2],dim(climD)[1],2)
p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[1,3],dim(climD)[1],2)
X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- exp(mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,]) + 
    sum(betas*Xtmp)
  tmp.mu <- tmp.mu + eta
  ex.mat[t,] <- rnbinom(ncol(ex.mat), mu=tmp.mu, size = mean_params[mean_params$Parameter=="phi","..1"])
}
rcp45 <- ex.mat


pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[2,2],dim(climD)[1],2)
p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[2,3],dim(climD)[1],2)
X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- exp(mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,]) + 
    sum(betas*Xtmp)
  tmp.mu <- tmp.mu + eta
  ex.mat[t,] <- rnbinom(ncol(ex.mat), mu=tmp.mu, size = mean_params[mean_params$Parameter=="phi","..1"])
}
rcp60 <- ex.mat


pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[3,2],dim(climD)[1],2)
p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[3,3],dim(climD)[1],2)
X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- exp(mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,]) + 
    sum(betas*Xtmp)
  tmp.mu <- tmp.mu + eta
  ex.mat[t,] <- rnbinom(ncol(ex.mat), mu=tmp.mu, size = mean_params[mean_params$Parameter=="phi","..1"])
}
rcp85 <- ex.mat


####
####  Plot climate change results
####
all.sims <- data.frame(pixel=c(1:length(eq.pixels)),
                       Baseline=eq.pixels,
                       RCP45=colMeans(rcp45),
                       RCP60=colMeans(rcp60),
                       RCP85=colMeans(rcp85))
all4plot <- melt(all.sims, id.vars = "pixel")

png("../results/climchange_small.png", width = 5, height=3, units = "in", res=200)
ggplot(all4plot)+
  geom_line(stat="density", aes(x=value, y=..density.., col=variable), size=1)+
  scale_color_manual(values=c("dodgerblue","tan","coral","darkred"), name="Scenario") +
  xlab("Equilibrium Percent Cover")+
  ylab("Estimated Kernel Density")+
  theme_bw()
dev.off()