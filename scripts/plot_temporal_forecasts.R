##  Script to plot results from all model temporal runs

##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 10-09-2015


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

longdf <- data.frame(model=NA,
                     scenario=NA,
                     param_set=NA,
                     timestep=NA,
                     cover=NA)

for(i in 1:length(all_files)){
  ttt <- readRDS(all_files[i])
  ttt2 <- adply(ttt, c(1,3))
  ttt3 <- melt(ttt2, id.vars = c("X1", "X2"))
  ttt3$variable <- as.numeric(ttt3$variable)
  colnames(ttt3) <- c("param_set", "id", "timestep", "cover")
  ttt3agg <- ddply(ttt3, .(param_set, timestep), summarise,
                   cover = mean(cover, rm.na=TRUE))
  ttt3agg$model <- strsplit(all_files[i], "_")[[1]][1]
  ttt3agg$scenario <- strsplit(all_files[i], "_")[[1]][2]
  tmpout <- ttt3agg[,c("model", "scenario", "param_set", "timestep", "cover")]
  longdf <- rbind(longdf, tmpout)
}
longdf <- longdf[2:nrow(longdf),]
longall <- longdf
longall$parammod <- paste(longall$param_set, longall$model)


# longdf <- subset(longall, model %in% c("bcc-csm1-1", "ccsm4"))
longdf <- longall
meandf <- ddply(longdf, .(timestep, scenario), summarise,
                avgcover=mean(cover, na.rm=TRUE),
                upcover=quantile(cover, 0.95, na.rm=TRUE),
                locover=quantile(cover, 0.05, na.rm=TRUE))
meandf$timestep <- rep(c(2011:2098), each=3)



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

mycols <- c("#F8DFC8",
            "#B2F395",
            "#F8BA61",
            "#68C9C8",
            "#A2AAD2",
            "#ED9991",
            "#AAB375",
            "#6CCA9C",
            "#FDCCF7",
            "#ECF7B5",
            "#D1F6DF",
            "#E0D979",
            "#97D1EC",
            "#F7B5C3",
            "#EDE2F1",
            "#A2C4A3",
            "#D3AE75",
            "#FAA873",
            "#D5AE9D",
            "#96BBBC",
            "#FEBB9B",
            "#B7AE55",
            "#BCADC0",
            "#EEFC8B",
            "#C0F0B3",
            "#E4AACF",
            "#BDAC8A",
            "#F2D088",
            "#87EFBB",
            "#ACF7E8",
            "#E5C963",
            "#CAAFB1",
            "#7AD5C2",
            "#D7BCE4",
            "#C0DCF0",
            "#ABBC8E",
            "#94BADD",
            "#FBE6E4",
            "#93D9E4",
            "#C6EE88",
            "#A2B6C5",
            "#E6A15F",
            "#E6F9EC",
            "#9DB870",
            "#C2B166",
            "#ACF4A4",
            "#F4C3A3",
            "#F3DAB6",
            "#DEFAD7",
            "#72DDB2",
            "#FBC867",
            "#E9DE89",
            "#E9AD91",
            "#F9F1F8",
            "#F3DCF0",
            "#74C995")
# ggplot()+
#   geom_line(data=obs_agg, aes(x=Year, y=avgcover))+
#   # geom_point(data=obs_agg, aes(x=Year, y=avgcover), size=4)+
#   geom_ribbon(data=meandf, aes(x=timestep, ymin=locover, ymax=upcover), alpha=0.7)+
#   # geom_line(data=longdf, aes(x=timestep, y=cover, color=model, group=parammod),alpha=0.5)+
#   facet_wrap("scenario")+
#   guides(color=FALSE)+
#   scale_y_continuous(limits=c(0,50))+
#   scale_color_manual(values=mycols)+
#   theme_bw()

ggplot()+
  geom_line(data=obs_agg, aes(x=Year, y=avgcover))+
  geom_ribbon(data=obs_agg, 
              aes(x=Year, ymin=locover, ymax=upcover),
              alpha=0.5)+
  geom_ribbon(data=meandf, 
              aes(x=timestep, ymin=locover, ymax=upcover, fill=scenario),
              alpha=0.5)+
  geom_line(data=meandf, aes(x=timestep, y=avgcover, color=scenario))+
  geom_vline(aes(xintercept=2011))+
  ylab("Sagebrush cover (%)")+
  guides(color=FALSE)+
  scale_y_continuous(limits=c(0,40))+
  scale_fill_manual(values=c("tan","coral","darkred"), name="IPCC \nScenario")+
  scale_color_manual(values=c("tan","coral","darkred"), name="IPCC \nScenario")+
  theme_bw()
ggsave("../temporal_forecast.png", width = 5, height = 4)


