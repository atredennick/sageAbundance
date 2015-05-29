##  Script to run GLMM fit on subset of sagebrush data

library(sageAbundance)
library(rstan)
library(ggmcmc)

####
##  Get data
####
datapath <-  "/Users/atredenn/Dropbox/sageAbundance_data/"
climD <- read.csv(paste(datapath,
                        "/studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv",
                        sep=""))
rawD <- read.csv(paste(datapath, "WY_SAGECoverData_V2check.csv", sep=""))

# Merge in climate data 
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Load knot data
load("../results/Knot_cell_distances.Rdata")


####
## Quick glm for initial values
####
mod <- glm(Cover ~ CoverLag, family="poisson", data=growD)
summary(mod)


####
##  Send data to stan function for fitting
####
mcmc <- model_nocovars(y = growD$Cover, lag = growD$CoverLag, K = K.data$K, 
                       cellid = growD$ID, iters = 50, warmup = 25, nchains = 1)



