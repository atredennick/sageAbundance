## Script to aggregate CMIP5 model projections for simulations

##  Relies on data downloaded from the CMIP5 website:
##  http://cmip-pcmdi.llnl.gov/cmip5/

##  Authors:      Andrew Tredennick and Peter Adler
##  Email:        atredenn@gmail.com
##  Date created: 10-10-2013



### Clean the workspace and set working dir
rm(list=ls()) 
setwd("../../data/climate/PROJECTIONS/CMIP5/")


####
####  Load necessary libraries -------------------------------------------------
####
library(reshape2)
library(grid)



####
####  Read in and formate precip projections -----------------------------------
####
ppt <- as.data.frame(read.csv("pr.csv", header=FALSE))
tmp<-read.table("COLS_pr.txt")
tmp<-as.character(tmp[,1])
colnames(ppt) <- c("Year", "Month",tmp)
ppt <- melt(ppt, id.vars=c("Year", "Month"))
ppt[,4] <- as.numeric(ppt[,4]) #NAs coerced for a couple December 2099 null values 
tmp<-unlist(strsplit(x=as.character(ppt$variable),split=".",fixed=T))
tmp<-matrix(tmp,nrow=length(tmp)/3,ncol=3,byrow=T)
colnames(tmp)<-c("model","rep","scenario")
ppt<-cbind(ppt,tmp)
ppt$period<-cut(ppt$Year,breaks=c(1950,2000,2050,2100),
  include.lowest=T,labels=c("past","present","future"))
ppt$season<-ifelse(ppt$Month > 6 & ppt$Month < 10,"summer","fall2spr")

pptMeans<-aggregate(value~Year+season+scenario+model,data=ppt,FUN=mean)

# Set number of days for each season
pptMeans$days<-ifelse(pptMeans$season=="summer",92,365-92)
# Multiply values (mean mm/day) by the number of days in season
pptMeans$ppt<-pptMeans$value*pptMeans$days




####
####  Read in and formate temperature projections ------------------------------
####
Tavg <- as.data.frame(read.csv("tas.csv", header=FALSE))
tmp<-read.table("COLS_tas.txt")
tmp<-as.character(tmp[,1])
colnames(Tavg) <- c("Year", "Month",tmp)
Tavg <- melt(Tavg, id.vars=c("Year", "Month"))
tmp<-unlist(strsplit(x=as.character(ppt$variable),split=".",fixed=T))
tmp<-matrix(tmp,nrow=length(tmp)/3,ncol=3,byrow=T)
colnames(tmp)<-c("model","rep","scenario")
Tavg<-cbind(Tavg,tmp)
Tavg$period<-cut(Tavg$Year,breaks=c(1950,2000,2050,2100),
  include.lowest=T,labels=c("past","present","future"))
Tavg$season<-ifelse(Tavg$Month > 3 & ppt$Month < 7,"spring","other")

TavgMeans<-aggregate(as.numeric(value)~Year+season+scenario+model,data=Tavg,FUN=mean)
colnames(TavgMeans) <- c("Year", "season", "scenario", "model", "TmeanSpr")



####
####  Merge and write files for use in simulations ---------------------------------------
####
tavg4out <- subset(TavgMeans, season!="other")
saveRDS(pptMeans, "../../../CMIP5_yearly_project_precipitation.RDS")
saveRDS(tavg4out, "../../../CMIP5_yearly_project_temperature.RDS")

