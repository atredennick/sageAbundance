##  Script to plot results from all model temporal runs

##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 10-09-2015


##  Clear the workspace...
rm(list=ls())


####
####  Load libraries -----------------------------------------------------------
####
library(ggplot2)
library(plyr)
library(reshape2)



####
####  Read in model simulations, 1 by 1, and stack -----------------------------
####
setwd("../results/yearlyforecasts/")
all_files <- list.files()
torms <- grep("MEAN", all_files)
all_files <- all_files[-torms]

longdf <- data.frame(model=NA,
                     scenario=NA,
                     paramset=NA,
                     year=NA,
                     cover=NA)

for(i in 1:length(all_files)){
  ttt <- readRDS(all_files[i])
  ttt2 <- as.data.frame(apply(ttt, c(2,1), mean))
  colnames(ttt2) <- paste0("params", 1:ncol(ttt2))
  ttt2$year <- c(2011:2098)
  ttt3 <- melt(ttt2, id.vars="year")
  ttt3$model <- strsplit(all_files[i], "_")[[1]][1]
  ttt3$scenario <- strsplit(all_files[i], "_")[[1]][2]
  colnames(ttt3) <- c("year", "paramset", "cover", "model", "scenario")
  tmpout <- ttt3[,c("model", "scenario", "paramset", "year", "cover")]
  longdf <- rbind(longdf, tmpout)
  print(paste(i, "out of", length(all_files)))
}


longdf <- longdf[2:nrow(longdf),]
longall <- longdf
longall$parammod <- paste(longall$param_set, longall$model)


# longdf <- subset(longall, model %in% c("bcc-csm1-1", "ccsm4"))
longdf <- longall
meandf <- ddply(longdf, .(year, scenario), summarise,
                avgcover=mean(cover, na.rm=TRUE),
                upcover=quantile(cover, 0.95, na.rm=TRUE),
                locover=quantile(cover, 0.05, na.rm=TRUE))



####
####  Read in a format observation data ----------------------------------------
####
obs_data <- read.csv("../../data/wy_sagecover_subset_noNA.csv")
obs_data <- subset(obs_data, Year>1984) # subsets out NA CoverLag values
obs_agg <- ddply(obs_data, .(Year), summarise,
                 avgcover=mean(Cover),
                 upcover=quantile(Cover, 0.95),
                 locover=quantile(Cover, 0.05))


####
####  Make the plot ------------------------------------------------------------
####
ggplot()+
  geom_line(data=obs_agg, aes(x=Year, y=avgcover))+
#   geom_ribbon(data=obs_agg, 
#               aes(x=Year, ymin=locover, ymax=upcover),
#               alpha=0.5)+
  geom_ribbon(data=meandf, 
              aes(x=year, ymin=locover, ymax=upcover, fill=scenario),
              alpha=0.5)+
  geom_line(data=meandf, aes(x=year, y=avgcover, color=scenario))+
  geom_vline(aes(xintercept=2011))+
  ylab("Mean sagebrush cover (%)")+
  guides(color=FALSE,fill=FALSE)+
  scale_y_continuous(limits=c(0,30))+
  scale_fill_manual(values=c("tan","coral","darkred"), 
                    name="IPCC \nScenario",
                    labels=c("RCP 4.5", "RCP 6.0", "RCP 8.5"))+
  scale_color_manual(values=c("tan","coral","darkred"), 
                     name="IPCC \nScenario",
                     labels=c("RCP 4.5", "RCP 6.0", "RCP 8.5"))+
  theme_bw()
ggsave("../temporal_forecast_wpois.png", width = 5, height = 4, dpi = 120)



####
####  PLOT A SHORT TERM FORECAST FROM ONE GCM
####
one_gcm <- subset(longdf, model == "miroc5" & year < 2022 & scenario=="rcp85")
ggplot()+
  geom_line(data=obs_agg, aes(x=Year, y=avgcover))+
  geom_line(data=one_gcm, aes(x=year, y=cover, group=paramset), alpha=0.5, color="darkred")+
  geom_vline(aes(xintercept=2011))+
  ylab("Mean sagebrush cover (%)")+
  guides(color=FALSE)+
  scale_y_continuous(limits=c(0,30))+
  theme_bw()
ggsave("../short_term_miroc.png", width = 5, height = 4, dpi=120)

