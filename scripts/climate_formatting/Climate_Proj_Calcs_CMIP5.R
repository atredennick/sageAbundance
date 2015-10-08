rm(list=ls()) #wipe workspace clean
#require(plyr)
require(reshape2)
require(grid)

## STUDY AREA 1 ##
setwd("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/climate/PROJECTIONS/CMIP5/")

## STUDY AREA 2 ##
# setwd("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea2/climate/PROJECTIONS/bcsd5/")


# Precip---------------------------------------------------------
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

# Temperature---------------------------------------------------------
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


pptMeans
TavgMeans

write.csv(pptMeans, "/Users/atredenn/Repos/sageAbundance/data/precipitation_projections.csv")
write.csv(TavgMeans, "/Users/atredenn/Repos/sageAbundance/data/temperature_projections.csv")


