##  Script to run GLMM fit on subset of sagebrush data

library(sageAbundance)
library(rstan)
library(ggmcmc)
library(parallel)

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
# growD <- subset(growD, Year<1994) # subset just 10 years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Load knot data
load("../results/Knot_cell_distances.Rdata")
K <- K.data$K

# Remove NA values
rms <- which(is.na(growD$Cover) == TRUE)
rmsK <- which(is.na(subset(growD, Year==1985)["Cover"])==TRUE)
K <- K[-rmsK,]
growD <- growD[-rms, ]
# if(length(which(is.na(growD$CoverLag)))>0){
#   rms <- which(is.na(growD$CoverLag) == TRUE)
#   growD <- growD[-rms, ]
#   K <- K[-rms,]
# }

# Set up new, continuous cell IDs
num_ids <- length(unique(growD$ID))
num_reps <- nrow(growD)/num_ids
new_ids <- rep(c(1:num_ids), num_reps)
growD$newID <- new_ids



####
## Quick glm for initial values
####
# mod <- glm(Cover ~ CoverLag, family="poisson", data=growD)
# summary(mod)


####
##  Send data to stan function for fitting
####
inits <- list()
inits[[1]] <- list(int_mu = 2, beta_mu = 0.05, 
                   alpha = rep(0,ncol(K)))
inits[[2]] <- list(int_mu = 2.5, beta_mu = 0.02, 
                   alpha = rep(0.5,ncol(K)))
inits[[3]] <- list(int_mu = 1.5, beta_mu = 0.09, 
                   alpha = rep(0.35,ncol(K)))

mcmc <- model_nocovars(y = growD$Cover, lag = growD$CoverLag, K = K, 
                       cellid = growD$newID, iters = 100, warmup = 50, 
                       nchains = 1, inits=inits)
saveRDS(mcmc, "stanmcmc_sage.RDS")


