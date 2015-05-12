# Bayesian model for sagebrush pixel growth

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

ptm <- proc.time()

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
rown <- 169
rows <- 1
ntest <- rown*rows

growD <- subset(fullD, ID<(ntest+1))
# growD <- subset(subD, is.na(CoverLag)==FALSE)
# growD <- subset(subD, Year>1984)
# tmpI <- which(growD$CoverLag==0)
# growD$CoverLag[tmpI] <- 0.001


####
#### Set up data structure for JAGS and iteration values
####
index <- which(is.na(growD$Cover)==FALSE) #vector of rows of non-missing data
N.obs <- growD$Cover[index]
lagObs <- growD$CoverLag[index]
cell  <- growD$ID[index]
nCells <- length(unique(cell))
cellMod <- c(1:nCells)
years <- growD$Year - 1983
nYears <- length(unique(years))
K <- K.data$K
nKnots <- ncol(K)

dataJ <- list(C = N.obs,
              lagC = lagObs,
              nObs = length(N.obs),
              cell = cell,
              nCells = nCells,
              cellMod = rep(cellMod, nYears),
              years = years,
              nYears = nYears,
              nKnots=nKnots,
              K=K,
              nTot = nrow(K))
n.adapt=100
n.update=200
n.iter=100
n.chains=1
n.thin=5

####
#### Run MCMC through JAGS
####
mod <- jags.model("sageGompertz_StateSpace_AST_JAGS.R", data=dataJ, n.chains=n.chains, n.adapt=n.adapt)
# update(mod, n.iter = n.update)
out <- coda.samples(mod, c("betaMu", "intMu", "alpha"),
                    n.iter=n.iter, n.thin=n.thin)
outStat <- summary(out)$stat  
write.csv(outStat, "outStats.csv")

print(proc.time() - ptm)


