#' Produce a climate dataframe for sageAbundance simulations
#' @param tdata subset of temperature dataframe for a single model and scenario
#' @param pdata subset of precipitation dataframe for a single model and scenario
#' @param years vector of desired years in output
#' @return clim_data formatted climate dataframe for simulating population model
#' @author Andrew Tredennick
format_climate <- function(tdata, pdata, years){
  ### Load 'plyr', if not already loaded
#   if("plyr" %in% rownames(installed.packages()) == FALSE){
#     stop("you need to install the 'plyr' package")
#   }
  if("package:plyr" %in% search() == FALSE) library(plyr)
  
  ### Add in 2 extra years to get lag precip != NA
  orig_years <- years                           # save original year vector
  years <- c(min(years)-2, min(years)-1, years) # add in 2 extra years
  
  ### Temperature
  t_tmp <- subset(tdata, Year %in% years)   # get appropriate years
  t_tmp <- t_tmp[ ,c("Year", "TmeanSpr")]   # subset relevant columns
  colnames(t_tmp) <- c("year", "TmeanSpr2") # rename TmeanSpr column
  t_tmp2 <- t_tmp                           # make duplicate dataframe
  t_tmp2$year <- t_tmp2$year+1              # add year for lag temperature
  colnames(t_tmp2)[2] <- "TmeanSpr1"        # rename for lag temperature
  t_tmp <- merge(t_tmp,t_tmp2,all.x=T)      # merge together
  
  ### Precipitation
  p_tmp <- subset(pdata, Year %in% years)      # get appropriate years
  p_tmp <- p_tmp[ ,c("Year", "season", "ppt")] # subset relevant columns
  
  # Lag (t-1) cumulative precipitation
  p_lag <- ddply(p_tmp, .(Year), summarise,
                 pptLag = sum(ppt))      # create summary df with cumulative ppt
  p_lag$Year <- p_lag$Year + 2           # make it lag by 2 years
  colnames(p_lag) <- c("year", "pptLag") # rename for consistency
  
  # Summer precipitation
  p_summer <- subset(p_tmp, season=="summer")     # subset for just summer
  p_summer <- p_summer[ ,c("Year", "ppt")]        # get relevant columns
  colnames(p_summer) <- c("year", "pptSummer2")   # rename ppt column
  p_summer2 <- p_summer                           # make duplicate dataframe
  p_summer2$year <- p_summer2$year + 1            # add year for lag precip
  colnames(p_summer2)[2] <- "pptSummer1"          # rename for lag precip
  p_summer <- merge(p_summer, p_summer2, all.x=T) # merge together
  
  # Spring precipitation
  p_spring <- subset(p_tmp, season=="fall2spr")   # subset for just spring
  p_spring <- p_spring[ ,c("Year", "ppt")]        # get relevant columns
  colnames(p_spring) <- c("year", "ppt2")         # rename ppt column
  p_spring2 <- p_spring                           # make duplicate dataframe
  p_spring2$year <- p_spring2$year + 1            # add year for lag precip
  colnames(p_spring2)[2] <- "ppt1"                # rename for lag precip
  p_spring <- merge(p_spring, p_spring2, all.x=T) # merge together
  
  p_tmp <- merge(p_spring, p_summer, all.x=T)     # merge both ppt dataframes
  p_tmp <- merge(p_tmp, p_lag, all.x=T)           # merge in pptLag
  
  ### Merge ppt and temp for output
  out_climate <- merge(t_tmp, p_tmp, all.x=T)            # merge all climate
  clim_data <- subset(out_climate, year %in% orig_years) # get original years
  return(clim_data)
}

## Function test
# all_temps <- readRDS("../data/CMIP5_yearly_project_temperature.RDS")
# all_ppts <- readRDS("../data/CMIP5_yearly_project_precipitation.RDS")
# 
# model_names <- unique(all_temps$model)
# tdata <- subset(all_temps, scenario=="rcp45" & model==model_names[1])
# pdata <- subset(all_ppts, scenario=="rcp45" & model==model_names[1])
# years <- c(2011:max(pdata$Year))
# 
# out_clim <- format_climate(tdata, pdata, years)
# head(out_clim)

