##  Script for leave-one-year-out cross validation of the
##  sagebrush abundance model. Here we use a simplified version
##  with no spatial random effect to reduce computational overhead.
##  We use a frequentist approach (lm).

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-09-2015

##  Clear the workspace...
rm(list=ls())



####
####  Load libraries -----------------------------------------------------------
####
library(ggplot2)
library(reshape2)
library(plyr)
library(lme4)



####
####  Read in observation data -------------------------------------------------
####
obs_data <- read.csv("../data/wy_sagecover_subset_noNA.csv")

# Get data structure right
growD <- subset(obs_data, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

# Get year information
all_years <- unique(growD$Year)
num_years <- length(all_years)

"%w/o%" <- function(x, y) x[!x %in% y] # x without y



####
####  Fit and predict total in sample observations -----------------------------
####
zids <- growD[which(growD$Cover==0),"ID"]
zids <- unique(c(zids, growD[which(growD$CoverLag==0),"ID"]))
pixels_to_keep <- growD$ID %w/o%  zids
modelD <- growD[which(growD$ID %in% pixels_to_keep),]

X <- modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X <- scale(X, center = TRUE, scale = TRUE)
modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")] <- X
# tmp.glm <- glmer(y~lag+X+(1|year), family="poisson")
tmp.lm <- glm(Cover~log(CoverLag)+pptLag+ppt1+ppt2+TmeanSpr1+TmeanSpr2, 
              family="poisson", data=modelD)
all.pred <- predict(tmp.lm, type="response")
sq.error.all <- as.numeric((all.pred-modelD$Cover)^2)
rmse.all <- sqrt(mean(sq.error.all))




####
####  Loop over leave out years and predict left out year ----------------------
####
rmse <- numeric(num_years)
counter <- 1
for(omit_year in all_years){
  tmpD <- subset(growD, Year != omit_year)
  # Remove pixels that have 0s in the time series
  zids <- tmpD[which(tmpD$Cover==0),"ID"]
  zids <- unique(c(zids, tmpD[which(tmpD$CoverLag==0),"ID"]))
  pixels_to_keep <- tmpD$ID %w/o%  zids
  modelD <- tmpD[which(tmpD$ID %in% pixels_to_keep),]
  obs_clim_means <- colMeans(modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")])
  obs_clim_sds <- apply(modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")], 2, sd)
  
  X <- modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  X <- scale(X, center = TRUE, scale = TRUE)
  modelD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")] <- X
  # tmp.glm <- glmer(y~lag+X+(1|year), family="poisson")
  tmp.lm <- glm(Cover~log(CoverLag)+pptLag+ppt1+ppt2+TmeanSpr1+TmeanSpr2, 
                family="poisson", data=modelD)
  
  ### Make predictions for omitted year
  predD <- subset(growD, Year == omit_year)
  climate_now <- predD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  # Scale climate predictors
  climate_now["pptLag"] <- (climate_now["pptLag"] - obs_clim_means["pptLag"])/obs_clim_sds["pptLag"]
  climate_now["ppt1"] <- (climate_now["ppt1"] - obs_clim_means["ppt1"])/obs_clim_sds["ppt1"]
  climate_now["ppt2"] <- (climate_now["ppt2"] - obs_clim_means["ppt2"])/obs_clim_sds["ppt2"]
  climate_now["TmeanSpr1"] <- (climate_now["TmeanSpr1"] - obs_clim_means["TmeanSpr1"])/obs_clim_sds["TmeanSpr1"]
  climate_now["TmeanSpr2"] <- (climate_now["TmeanSpr2"] - obs_clim_means["TmeanSpr2"])/obs_clim_sds["TmeanSpr2"]
  predD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")] <- climate_now
  tmp.pred <- predict(tmp.lm, newdata = predD, type="response")
  sq.error <- as.numeric((tmp.pred-predD$Cover)^2)
  rmse[counter] <- sqrt(mean(sq.error))
  print(omit_year)
  counter <- counter+1
}
mean(rmse)
hist(rmse)

