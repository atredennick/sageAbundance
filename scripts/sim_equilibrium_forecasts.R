##
##  R Script to Simulate Equilibrium Forecasts Under Projected Climate Change
##
##  Climate data comes from the CMIP5 projections
##  We simulate equilibrium cover for projected changes by the 2050-2100
##  period for each climate model. Then we average the results over climate
##  models for each RCP scenario.
##
##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 12-7-2015
##

# Clear the Workspace
rm(list=ls())



####
####  Set File Paths
####
datapath <- "../data/"
knotpath <- "../results/"



####
####  Load Libraries
####
library(rstan)
library(reshape2)
library(ggplot2)
library(plyr)
library(sageAbundance)
library(gridExtra)
library(RColorBrewer)



####
####  Set Plotting Themes
####
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=16),text=element_text(size=16),
                legend.text=element_text(size=16))
tmp.theme2=theme(strip.text=element_text(face="bold"),
                 axis.title=element_text(size=16),text=element_text(size=16),
                 legend.text=element_text(size=16))
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))



####
####  Load Observed Data, Climate Data, Climate Projections, and MCMC Results
####
fullD <- read.csv(paste0(datapath,"wy_sagecover_subset_noNA.csv"))
if(length(which(is.na(fullD$Cover))) > 0) stop("data contains NA values")

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Climate data
climD <- read.csv(paste0(datapath,"/climate/DAYMET/FormattedClimate_WY_SA1.csv"))

# Knots and MCMC results
load(paste0(knotpath,"Knot_cell_distances_subset.Rdata"))
K <- K.data$K
outs <- readRDS("../results/poissonSage_randYear_mcmc.RDS")



####
####  Fit Colonization Logistic Model for 0% Cover Cells
####
colD <- growD[which(growD$CoverLag == 0), ]
colD$colonizes <- ifelse(colD$Cover==0, 0, 1)
col.fit <- glm(colonizes ~ 1, data=colD, family = "binomial")
col.intercept <- as.numeric(coef(col.fit))
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
avg.new.cover <- round(mean(colD[which(colD$Cover>0),"Cover"]),0)



####
####  Define Simulation Function
####
sim_sage <- function(init_cover_vector, scaled_clim_matrix, 
                     avg_new_cover, params, time_steps){
  int_mu <- params[["int_mu"]]
  beta_mu <- params[["beta_mu"]]
  betas <- params[["betas"]]
  eta <- params[["eta"]]
  cover_matrix <- matrix(NA,nrow=time_steps,ncol=length(init_cover_vector))
  cover_matrix[1,] <- init_cover_vector
  X <- scaled_climate_matrix
  for(t in 2:time_steps){
    ## Growth Process
    Xtmp <- X[sample(c(1:nrow(X)), 1),]
    dens.dep <- beta_mu*log(cover_matrix[t-1,])
    dens.dep[which(dens.dep==-Inf)] <- 0
    tmp.mu <- int_mu + dens.dep + sum(betas*Xtmp)
    tmp.mu <- exp(tmp.mu + eta)
    tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
    
    ## Colonization Process
    zeros <- which(ex.mat[t-1,]==0)
    colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
    colonizer.cover <- colonizers*avg_new_cover
    tmp.out[zeros] <- colonizer.cover
    ex.mat[t,] <- tmp.out
  } # end simulation loop
  return(ex.mat)
} # end simulation function



####
####  Run equilibrium population simulation ------------------------------------
####
time.steps <- 2000
burn.in <- 100
mean_params <- ddply(outs, .(Parameter), summarise,
                     value = mean(value))
alphas <- mean_params[grep("alpha", mean_params$Parameter),"value"]
betas <- mean_params[grep("beta", mean_params$Parameter),"value"][2:6]
eta <- K%*%alphas
pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
clim_sim <- climD[climD$year %in% unique(growD$Year),]
X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)

for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  dens.dep <- mean_params[mean_params$Parameter=="beta_mu","value"]*log(ex.mat[t-1,])
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","value"] + dens.dep + sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  
  #Colonization
  zeros <- which(ex.mat[t-1,]==0)
  colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
  colonizer.cover <- colonizers*avg.new.cover
  tmp.out[zeros] <- colonizer.cover
  
  ex.mat[t,] <- tmp.out
}

