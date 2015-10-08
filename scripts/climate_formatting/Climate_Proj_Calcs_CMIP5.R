## Script to average and aggregate all CMIP5 model projections

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
ppt <- melt(ppt, id.vars=c("Year", "Month"),)
ppt[,4] <- as.numeric(ppt[,4]) #NAs coerced for a couple December 2099 null values 
tmp<-unlist(strsplit(x=as.character(ppt$variable),split=".",fixed=T))
tmp<-matrix(tmp,nrow=length(tmp)/3,ncol=3,byrow=T)
colnames(tmp)<-c("model","rep","scenario")
ppt<-cbind(ppt,tmp)
ppt$period<-cut(ppt$Year,breaks=c(1950,2000,2050,2100),
  include.lowest=T,labels=c("past","present","future"))
ppt$season<-ifelse(ppt$Month > 6 & ppt$Month < 10,"summer","fall2spr")

pptMeans<-aggregate(value~period+season+scenario,data=ppt,FUN=mean)
colnames(pptMeans) <- c("period","season","scenario","value")
pptMeans$days<-ifelse(pptMeans$season=="summer",92,365-92)
pptMeans$ppt<-pptMeans$value*pptMeans$days
pptMeans<-reshape(pptMeans[,c("period","season","scenario","ppt")],
              idvar=c("season","scenario"),timevar="period",direction="wide")
pptMeans$change<-(pptMeans$ppt.future-pptMeans$ppt.past)/pptMeans$ppt.past

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

TavgMeans<-aggregate(as.numeric(value)~period+season+scenario,data=Tavg,FUN=mean)
colnames(TavgMeans) <- c("period","season","scenario","value")
TavgMeans<-reshape(TavgMeans[,c("period","season","scenario","value")],
              idvar=c("season","scenario"),timevar="period",direction="wide")
TavgMeans$change<-TavgMeans$value.future-TavgMeans$value.past



####
####  Write files for use in simulations ---------------------------------------
####
write.csv(pptMeans, "../../../precipitation_projections.csv")
write.csv(TavgMeans, "../../../temperature_projections.csv")


