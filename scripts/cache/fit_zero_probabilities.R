##  Script to run GLM fit on colonization/extinction probabilities
##  Uses simple lm call
##
##  Author:      Andrew Tredennick
##  Email:       atredenn@gmail.com
##  Last Update: 9-29-2015


####
####  Set some file paths, etc.
####
datapath <- "/Users/atredenn/Dropbox/sageAbundance_data/"



####
####  Get data
####
fullD <- read.csv(paste0(datapath,"wy_sagecover_subset_noNA.csv"))
if(length(which(is.na(fullD$Cover))) > 0) stop("data contains NA values")

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data



####
####  Fit colonization logistic ------------------------------------------------
####
colD <- growD[which(growD$CoverLag == 0), ]
colD$colonizes <- ifelse(colD$Cover==0, 0, 1)

col.fit <- glm(colonizes ~ 1, data=colD, family = "binomial")
col.intercept <- as.numeric(coef(col.fit))


####
####  Test the fit for probability of colonization -----------------------------
####
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
# rbinom(100,1,antilogit(col.intercept))
avg.new.cover <- round(mean(colD[which(colD$Cover>0),"Cover"]),0)
