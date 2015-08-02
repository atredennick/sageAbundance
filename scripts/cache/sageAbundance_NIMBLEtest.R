#clear everything, just to be safe 
rm(list=ls(all=TRUE))

####
#### Load libraries
####
library(nimble)
library(ggplot2)
library(tidyr)

####
#### Bring in data
####
studyArea="1"
ifelse(studyArea=="1",
       climD<-read.csv("../../studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv"),
       climD<-read.csv("../../studyarea2/climate/DAYMET/FormattedClimate_WY_SA2.csv"))
ifelse(studyArea=="1",
       rawD<-read.csv("../../cover_structure/WY_SAGECoverData_V2check.csv"),
       rawD<-read.csv("../../cover_structure/WY_SAGECoverData_SA2.csv"))

# merge in climate data 
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)

# get knot data
load("Knot_cell_distances.Rdata")

####
#### Take subset of data for test fitting
####
# rown <- 169
# rows <- 1
# ntest <- rown*rows
# 
# subD <- subset(fullD, ID<(ntest+1))
# growD <- subset(subD, is.na(CoverLag)==FALSE)
growD <- subset(fullD, Year>1984)
growD <- subset(growD, Year<1986)


# growD <- subset(fullD, Year>=1985 & Year<=1987)
# coords.mod <- as.matrix(growD[,c("Lon", "Lat")])
# D <- as.matrix(dist(coords.mod))

####
#### Quick glm for initial values ----------------------
####
mod <- glm(Cover ~ CoverLag, family="poisson", data=growD)
summary(mod)
init_alpha <- as.numeric(coef(mod)[1])
init_beta <- as.numeric(coef(mod)[2])

####
#### Set up data structure for NIMBLE ----------------------
####
index <- which(is.na(growD$Cover)==FALSE) #vector of rows of non-missing data
N.obs <- growD$Cover[index]
lagObs <- growD$CoverLag[index]
cell  <- growD$ID[index]
K <- K.data$K[index,]
nKnots <- ncol(K)

constants <- list(nObs = length(N.obs), 
                  nKnots = nKnots,
                  N = lagObs,
                  K = K,
                  cell = cell)
constants2 <- list(nObs = length(N.obs), 
                  N = lagObs)
data <- list(C = N.obs)
inits <- list(intMu=rnorm(1,init_alpha,0.05), 
              betaMu=rnorm(1,init_beta,0.05), 
              alpha=rep(0,nKnots))

####
#### Set up the BUGS model ----------------------------
####
ast <- nimbleCode({
  for(i in 1:nObs){
    C[i] ~ dpois(lambda[i])
    lambda[i] <- exp(intMu + betaMu*N[i] + eta[cell[i]]) 
  }
  eta[1:nObs] <- K[,]%*%alpha[]
  for(j in 1:nKnots){
    alpha[j] ~ dnorm(0,tau)
  }
  intMu ~ dunif(-5,5)
  betaMu ~ dunif(0, 1)
  tau ~ dgamma(1,0.01)
})

# ast <- nimbleCode({
#   for(i in 1:nObs){
#     C[i] ~ dpois(lambda[i])
#     lambda[i] <- exp(intMu + betaMu*N[i])
#   }
#   intMu ~ dnorm(1.23,1.9)
#   betaMu ~ dnorm(0.08, 0.1)
# })

####
#### Compile and run NIMBLE sampler ---------------------
####
ptm <- proc.time()
Rmodel <- nimbleModel(code = ast, 
                      constants = constants, 
                      data = data,
                      dimensions = list(K = c(dim(K)[1], dim(K)[2]),
                                        alpha = c(dim(K)[2])))

mcmcspec <- configureMCMC(Rmodel, print=TRUE, thin=100)
mcmcspec$addSampler('RW_block', list(targetNodes = c('intMu', 'betaMu'),
                                     adaptInterval = 100))
mcmcspec$addMonitors(c('alpha', 'betaMu', 'intMu', 'tau'))
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

Cmcmc$run(100000)
# Cmcmc$run(5000)

####
#### Collect, save, and plot samples ------------------------
####
samples <- as.data.frame(as.matrix(Cmcmc$mvSamples))
saveRDS(samples, "MCMCsamplesNIMBLE.rds")
proc.time() - ptm

# samples <- readRDS("MCMCsamplesNIMBLE.rds")
# dim(samples)
# samples1 <- samples[,c("betaMu", "intMu", "tau", "alpha[12]")]
# long <- gather(samples1)
# apply(samples1, 2, mean)
# 
# ggplot(long) + 
#   geom_line(aes(seq_along(value), value)) + 
#   facet_wrap(~key, scale='free')
# 
# ggplot(long) + 
#   geom_density(aes(value)) + 
#   facet_wrap(~key, scale='free')



