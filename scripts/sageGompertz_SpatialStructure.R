# Bayesian model for sagebrush pixel growth

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)

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
# rows <- 40
# ntest <- rown*rows
# 
# subD <- subset(fullD, ID<(ntest+1))
# growD <- subset(subD, is.na(CoverLag)==FALSE)
growD <- subset(fullD, Year>1984)
# tmpI <- which(growD$CoverLag==0)
# growD$CoverLag[tmpI] <- 0.001

# growD <- subset(fullD, Year>=1985 & Year<=1987)
# coords.mod <- as.matrix(growD[,c("Lon", "Lat")])
# D <- as.matrix(dist(coords.mod))

####
#### Set up data structure for JAGS and iteration values
####
index <- which(is.na(growD$Cover)==FALSE) #vector of rows of non-missing data
N.obs <- growD$Cover[index]
lagObs <- growD$CoverLag[index]
cell  <- growD$ID[index]
K <- K.data$K
nKnots <- ncol(K)

dataJ <- list(C = N.obs,
              N = lagObs,
              nObs = length(N.obs),
              cell = cell,
              nKnots=nKnots,
              K=K,
              nTot = nrow(K))
n.adapt=50
n.update=1
n.iter=10
n.chains=1
n.thin=0

####
#### Run MCMC through JAGS
####
ptm <- proc.time()
mod <- jags.model("sageGompertz_SpatialStructure_JAGS.R", data=dataJ, n.chains=n.chains, n.adapt=n.adapt)
# update(mod, n.iter = n.update)
out <- coda.samples(mod, c("betaMu", "intMu", "alpha"),
                    n.iter=n.iter, n.thin=n.thin)
print(proc.time() - ptm)

outStat <- summary(out)$stat  
write.csv(outStat, "outStats.csv")


