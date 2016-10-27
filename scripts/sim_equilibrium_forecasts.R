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

global_time_steps <- 500
global_burnin <- 100


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

# Climate projections
ppt_projs <- subset(read.csv("../data/precipitation_projections_bymodel.csv"), 
                    scenario!="rcp26" & season=="fall2spr")
temp_projs <- subset(read.csv("../data/temperature_projections_bymodel.csv"),
                     scenario!="rcp26" & season=="spring")
projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
                  "deltaPpt"=1+ppt_projs$change,
                  "deltaTspr"=temp_projs$change,
                  "model"=ppt_projs$model)

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
  int_yrs <- params[["int_yrs"]]
  beta_mu <- params[["beta_mu"]]
  betas <- params[["betas"]]
  eta <- params[["eta"]]
  cover_matrix <- matrix(NA,nrow=time_steps,ncol=length(init_cover_vector))
  cover_matrix[1,] <- init_cover_vector
  X <- scaled_climate_matrix
  pb <- txtProgressBar(min=2, max=time_steps, char="+", style=3, width=65)
  for(t in 2:time_steps){
    ## Growth Process
    int_now <- int_yrs[sample(c(1:length(int_yrs)), 1)]
    Xtmp <- X[sample(c(1:nrow(X)), 1),]
    dens.dep <- beta_mu*log(cover_matrix[t-1,])
    dens.dep[which(dens.dep==-Inf)] <- 0
    tmp.mu <- int_now + dens.dep + sum(betas*Xtmp)
    tmp.mu <- exp(tmp.mu + eta)
    tmp.out <- rpois(ncol(cover_matrix), lambda = tmp.mu)
    
    ## Colonization Process
    zeros <- which(cover_matrix[t-1,]==0)
    colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
    colonizer.cover <- colonizers*avg_new_cover
    tmp.out[zeros] <- colonizer.cover
    cover_matrix[t,] <- tmp.out
    setTxtProgressBar(pb, t)
  } # end simulation loop
  return(cover_matrix)
} # end simulation function



####
####  Set Up Parameter List for Simulations
####
mean_params <- ddply(outs, .(Parameter), summarise,
                     value = mean(value))
alphas <- mean_params[grep("alpha", mean_params$Parameter),"value"]
betas <- mean_params[grep("beta", mean_params$Parameter),"value"][2:6]
eta <- K%*%alphas
beta_mu <- mean_params[mean_params$Parameter=="beta_mu","value"]
# int_mu <- mean_params[mean_params$Parameter=="int_mu","value"]
yrint_ids <- grep("int_yr", mean_params$Parameter)
int_yrs <- mean_params[yrint_ids,"value"]
params <- list(int_yrs=int_yrs, beta_mu=beta_mu, betas=betas, eta=eta)



####
####  Run Equilibrium Population Simulation
####
burn_in <- global_burnin
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
init_cover_vector <- rep(1, times=nrow(subset(growD, Year==1985)))
Xtmp <- climD[climD$year %in% unique(growD$Year),clim_vars]
scaled_climate_matrix <- scale(Xtmp, center = TRUE, scale = TRUE)
avg_new_cover <- avg.new.cover
time_steps <- global_time_steps
equil_simulation <- sim_sage(init_cover_vector, scaled_clim_matrix, 
                             avg_new_cover, params, time_steps)
equil_results <- colMeans(equil_simulation[(burn_in+1):time_steps,])



####
####  Run Climate Change Simulations by Model
####
obs_clim<-climD[climD$year %in% unique(growD$Year),clim_vars]
clim_avg <- apply(X = obs_clim, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = obs_clim, MARGIN = 2, FUN = sd)

sim_matrix <- matrix(NA, nrow = nrow(projC), ncol = nrow(subset(growD, Year==1985)))
time_steps <- global_time_steps
burn_in <- global_burnin
for(i in 1:nrow(projC)){
  p.climD <- obs_clim
  p.climD[,c(1:3)]<-obs_clim[,c(1:3)]*matrix(projC[i,"deltaPpt"],dim(climD)[1],3)
  p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[i,"deltaTspr"],dim(climD)[1],2)
  X_sim = p.climD[,clim_vars]
  
  # Now scale based on perturbed or regular data, depending on scenario
  X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
  X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
  X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
  X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
  X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
  scaled_climate_matrix <- X_sim
  current_simulation <- sim_sage(init_cover_vector, scaled_clim_matrix, 
                               avg_new_cover, params, time_steps)
  current_results <- colMeans(current_simulation[(burn_in+1):time_steps,])
  sim_matrix[i,] <- current_results
  print(i)
}



####
####  Average Simulation Results by RCP Scenario
####
rcps <- c("rcp45", "rcp60", "rcp85")
rcp_matrix <- matrix(NA, nrow = length(rcps), ncol=ncol(sim_matrix))
for(i in 1:length(rcps)){
  ids <- which(projC$scenario==rcps[i])
  rcp_matrix[i,] <- colMeans(sim_matrix[ids,])
}



####
####  Combine Output and Make Spatil Plot
####
proj.equil <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
                         Lat=subset(growD, Year==1985)$Lat,
                         CURRENT=equil_results,
                         RCP45=rcp_matrix[1,],
                         RCP60=rcp_matrix[2,],
                         RCP85=rcp_matrix[3,])
colnames(proj.equil)[4:6] <- c("RCP 4.5", "RCP 6.0", "RCP 8.5")
proj.equil2 <- melt(proj.equil, id.vars = c("Lon", "Lat"))

# png("../results/climchange_small.png", width = 6, height=4.5, units = "in", res=150)
ggplot(proj.equil2, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=value, fill=value))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover")+
  facet_wrap("variable", ncol=2)+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))
ggsave("../results/clim_change_mean_spatial.png", height=8, width=8)





