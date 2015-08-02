##  Script to run GLM fit on subset of sagebrush data to test
##  effect of log vs. non-log lag cover

# Clear workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(sageAbundance)
library(rstan)
library(ggmcmc)
library(parallel)

####
####  Get data
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


####
####  Test log transform effect in GLM
####
y = growD$Cover
lag = growD$CoverLag
K = K.data$K
cellid = growD$ID
X = growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X = scale(X, center = TRUE, scale = TRUE)

mod <- glm(y ~ (lag+1)+X, family="poisson", data=growD)
mod.log <- glm(y ~ log(lag+1)+X, family="poisson", data=growD)
summary(mod)
summary(mod.log)
library(MASS)
summary(m1 <- glm.nb(y ~ (lag+1)+X, data = growD))
pchisq(2 * (logLik(m1) - logLik(mod)), df = 1, lower.tail = FALSE)

####
####  Simulate the population (poisson)
####
time.steps <- 100
pixels <- 100
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
clim_sim <- climD[climD$year %in% unique(growD$Year),]
X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- exp(coef(mod)[1] + coef(mod)[2]*ex.mat[t-1,]) + sum(coef(mod)[3:7]*Xtmp)
  ex.mat[t,] <- rpois(ncol(ex.mat), tmp.mu)
}

matplot(c(1:time.steps), ex.mat, type="l", col="grey")
hist(y, freq = FALSE)
lines(density(ex.mat, adjust = 4), col="red", lwd=2)
mean(y)
mean(ex.mat)



####
####  Simulate the population (negative binomial)
####
time.steps <- 1000
pixels <- 100
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
clim_sim <- climD[climD$year %in% unique(growD$Year),]
X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- exp(coef(m1)[1] + coef(m1)[2]*ex.mat[t-1,]) + sum(coef(m1)[3:7]*Xtmp)
  ex.mat[t,] <- rnbinom(ncol(ex.mat), mu=tmp.mu, size = m1$theta)
}

matplot(c(1:time.steps), ex.mat, type="l", col="grey")
hist(y, freq = FALSE)
lines(density(ex.mat, adjust = 4), col="red", lwd=2)
mean(y)
mean(ex.mat)
median(y)
median(ex.mat)




