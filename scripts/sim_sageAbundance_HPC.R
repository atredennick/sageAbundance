###########################################################
##  Supporting Material for Tredennick et al. 201x       ##  
##  Code to reproduce results displayed in figures x,x,x ##
###########################################################

##  Script to simulate fitted sage population abundance model
##  using observed and projected climate. This script runs simulations
##  for each MCMC iteration combination of parameters. That way we get
##  the full uncertainty around model forecasts.

##  Since this is a large job, this script is run on the Utah State
##  University super computer. Thus, file paths are different from
##  other scripts. This also means that you should probably not run
##  this as is on your PC. It will take forever...

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-06-2015



####
#### Clear workspace -----------------------------------------------------------
####
rm(list=ls(all=TRUE))



####
####  Load required libraries --------------------------------------------------
####
library(rstan)
library(reshape2)
library(plyr)


####
####  Set global simulation parameters -----------------------------------------
####
time.steps <- 600



####
####  Set climate driver id programmatically -----------------------------------
####
##  For HPC runs
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
climate_id <- as.numeric(myargument)

##  For PC tests (comment out for HPC runs!!)
# climate_id <- 1



####
####  Set some file paths, etc. ------------------------------------------------
####
##  For HPC runs
datapath <- ""
datapath2 <- ""
resultspath <- ""

##  For PC runs (comment out for HPC runs!!)
# datapath <- "/Users/atredenn/Dropbox/sageAbundance_data/"
# datapath2 <- "../data/"
# resultspath <- "../results/"



####
####  Get data -----------------------------------------------------------------
####
# Read in full data set
fullD <- read.csv(paste0(datapath,"wy_sagecover_subset_noNA.csv"))
# Make sure there aren't any leftover NAs, or wrong data set
if(length(which(is.na(fullD$Cover))) > 0){
  stop("data contains NA values")
}

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Get observed climate time series
climD <- read.csv(paste0(datapath2,"FormattedClimate_WY_SA1.csv"))

# Get projected climate changes
ppt_projs <- read.csv(paste0(datapath2, "precipitation_projections.csv"))
ppt_projs <- subset(ppt_projs, season=="fall2spr" & scenario != "rcp26")
temp_projs <- read.csv(paste0(datapath2, "temperature_projections.csv"))
temp_projs <- subset(temp_projs, season=="spring" & scenario != "rcp26")
scenarios <- as.character(ppt_projs$scenario)
climate_projections <- data.frame("scenario"= c("obs", scenarios),
                                  "deltaPpt" = c(1,ppt_projs$change+1),
                                  "deltaTspr"= c(0,temp_projs$change))



####
####  Load knot data and MCMC results ------------------------------------------
####
load(paste0(resultspath,"Knot_cell_distances_subset.Rdata"))
K <- K.data$K
outs <- readRDS(paste0(resultspath,"poissonSage_randYear_mcmc.RDS"))
colnames(outs)[which(colnames(outs)=="chain")] <- "Chain" # for clarity in loop


####
####  Fit colonization logistic ------------------------------------------------
####
colD <- growD[which(growD$CoverLag == 0), ]
colD$colonizes <- ifelse(colD$Cover==0, 0, 1)
col.fit <- glm(colonizes ~ 1, data=colD, family = "binomial")
col.intercept <- as.numeric(coef(col.fit))
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
avg.new.cover <- round(mean(colD[which(colD$Cover>0),"Cover"]),0)



####
####  Set up scaled climate for simulation -------------------------------------
####
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
p.climD <- climD[climD$year %in% unique(growD$Year), clim_vars]

# Subset climate scenario
clim_change <- climate_projections[climate_id,]
change_matrix_ppt <- matrix(as.numeric(clim_change["deltaPpt"]),
                            dim(p.climD)[1],2)
change_matrix_temp <- matrix(as.numeric(clim_change["deltaTspr"]),
                             dim(p.climD)[1],2)

# Climate perturbations
clim_avg <- apply(X = p.climD, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = p.climD, MARGIN = 2, FUN = sd)
p.climD[,c(2:3)] <- p.climD[,c(2:3)] * change_matrix_ppt
p.climD[,c(4:5)] <- p.climD[,c(4:5)] + change_matrix_temp
X_sim = p.climD[,clim_vars]

# Now scale based on perturbed or regular data, depending on scenario
X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]



####
####  Set up model and other parameters for simulation -------------------------
####
# Extract spatial field effects for eta calculation
alpha_data <- outs[grep("alpha", outs$Parameter),]

# Make vector for extracting climate effect names in MCMC dataframe
climeffs <- c("beta.1.", "beta.2.", "beta.3.", "beta.4.", "beta.5.")

# Set number of chains and iterations
nchains <- length(unique(outs$Chain))
niters <- length(unique(outs$mcmc_iter))
totsims <- nchains*niters

# Get number of pixels to simulate
pixels <- nrow(subset(growD, Year==1985))

# Make empty array for full set of simulations
ex.arr <- array(NA, dim=c(totsims, time.steps, pixels))
ex.arr[,1,] <- 1 # initialize all pixels at 1% cover



####
####  Main simulation loop -----------------------------------------------------
####
counter <- 1 # for keeping track of full number of iterations
for(i in 1:nchains){ # loop through chains
  
  for(j in 1:niters){ # loop through iterations within chain
    
    chain <- i
    iter <- j
    int_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="int_mu")[,"value"])
    beta_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="beta_mu")[,"value"])
    betas <- subset(outs, Chain==chain & mcmc_iter==iter & Parameter%in%climeffs)[,"value"]
    betas <- as.numeric(unlist(betas))
    alphas <- subset(alpha_data, Chain==chain & mcmc_iter==iter)[,"value"]
    eta <- K%*%as.numeric(unlist(alphas))
    
    for(t in 2:time.steps){ # loop throug time steps for current parameter set
      Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
      dens.dep <- beta_mu*log(ex.arr[counter,t-1,])
      tmp.mu <- int_mu + dens.dep + sum(betas*Xtmp)
      tmp.mu <- exp(tmp.mu + eta)
      tmp.out <- rpois(pixels, lambda = tmp.mu)
      
      #Colonization
      zeros <- which(ex.arr[counter,t-1,]==0)
      colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
      colonizer.cover <- colonizers*avg.new.cover
      tmp.out[zeros] <- colonizer.cover
      
      ex.arr[counter,t,] <- tmp.out
    } # next time step
    
  counter <- counter+1  
  
  } # next MCMC iteration
  
} # next MCMC chain



####
####  Format and save output ---------------------------------------------------
####
# Make the output a list: dim[1] = parameter set
#                         dim[2] = time step
#                         dim[3] = pixel
output <- alply(ex.arr, 1) 
filename <- paste0(climate_projections[climate_id,"scenario"], 
                   "clim_paramvary_sims.RDS")
saveRDS(output, filename)



