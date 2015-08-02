####
####  Read in Data
####

library(raster)
climD <- read.csv("FormattedClimate_WY_SA1.csv")
rawD<-read.csv("WY_SAGECoverData_V2check.csv")
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)
head(fullD)

years=sort(unique(fullD$Year))
n.years=length(years)
sites=sort(unique(fullD$ID))
n.sites=length(sites)

gray.colors.rev <- function (n, start = 0, end = 1, gamma = 1){ 
  gray(seq.int(to= start^gamma, from = end^gamma, length.out = n)^(1/gamma))
}

####
####  Make Initial Maps 
####

t=1
year.tmp=(fullD$Year==years[t])
cover.tmp=fullD[year.tmp,'Cover']
lon.locs=sort(unique(fullD[year.tmp,'Lon']))
lat.locs=sort(unique(fullD[year.tmp,'Lat']))
n.lon=length(lon.locs)
n.lat=length(lat.locs)
#plot(fullD[year.tmp,c('Lon','Lat')],cex=.2)
image(matrix(cover.tmp,n.lon,n.lat),x=lon.locs,y=lat.locs,asp=TRUE)
points(sample(x = lon.locs, size = 400), sample(lat.locs, 400), pch=19)





cover.mean.vec=rep(0,n.years)
cover.mat=matrix(0,n.sites,n.years)
png(file="sage_maps.png",width=10,height=5,units="in",res=200)
layout(matrix(1:n.years,4,7,byrow=TRUE))
par(xaxt="n",yaxt="n",mar=c(1.5,1,2,1))
for(t in 1:n.years){
  year.tmp=(fullD$Year==years[t])
  cover.tmp=fullD[year.tmp,'Cover']
  cover.mat[,t]=cover.tmp
  cover.mean.vec[t]=mean(cover.tmp,na.rm=TRUE)
  lon.locs=sort(unique(fullD[year.tmp,'Lon']))
  lat.locs=sort(unique(fullD[year.tmp,'Lat']))
  n.lon=length(lon.locs)
  n.lat=length(lat.locs)
  image(matrix(cover.tmp,n.lon,n.lat),x=lon.locs,y=lat.locs,asp=TRUE,col=gray.colors.rev(100),main=years[t],zlim = c(0,35))
}
dev.off()


####
#### Variograms for each year
####
library(gstat)



####
####  Check for homogeneous shift in value
####

plot(years,cover.mean.vec,type="o")


####
####  Plot Difference Maps 
####

diff.mat=matrix(0,n.sites,n.years)
png(file="sage_diff_maps.png",width=10,height=5,units="in",res=200)
layout(matrix(1:n.years,4,7,byrow=TRUE))
par(xaxt="n",yaxt="n",mar=c(1.5,1,2,1))
for(t in 2:n.years){
  year.1.tmp=(fullD$Year==years[t-1])
  year.2.tmp=(fullD$Year==years[t])
  cover.1.tmp=fullD[year.1.tmp,'Cover']
  cover.2.tmp=fullD[year.2.tmp,'Cover']
  cover.diff.tmp=cover.2.tmp-cover.1.tmp
  diff.mat[,t]=cover.diff.tmp
  lon.locs=sort(unique(fullD[year.tmp,'Lon']))
  lat.locs=sort(unique(fullD[year.tmp,'Lat']))
  n.lon=length(lon.locs)
  n.lat=length(lat.locs)
  image(matrix(cover.diff.tmp,n.lon,n.lat),x=lon.locs,y=lat.locs,col=gray.colors.rev(100),asp=TRUE,main=c(years[t-1],years[t]))
}
dev.off()

matplot(diff.mat[,13:15],type="l",lty=1)
abline(h=0,col=8)

layout(matrix(1:3,3,1))
plot(diff.mat[,13],type="l",lty=1,col=1)
abline(h=0,col=8)
plot(diff.mat[,14],type="l",lty=1,col=2)
abline(h=0,col=8)
plot(diff.mat[,15],type="l",lty=1,col=3)
abline(h=0,col=8)

ylim=c(-20,20)
xlim=c(0,40)
layout(matrix(1:3,1,3))
plot(cover.mat[,12],diff.mat[,13],col=1,ylim=ylim,xlim=xlim)
abline(0,-1)
plot(cover.mat[,13],diff.mat[,14],col=2,ylim=ylim,xlim=xlim)
abline(0,-1)
plot(cover.mat[,14],diff.mat[,15],col=3,ylim=ylim,xlim=xlim)
abline(0,-1)

plot(years,apply(diff.mat,2,mean,na.rm=TRUE),type="o")


####
####  Center Each Year 
####

cover.mean.vec=rep(0,n.years)
cover.mean.center.vec=rep(0,n.years)
png(file="sage_center_maps.png",width=10,height=5,units="in",res=200)
layout(matrix(1:n.years,4,7,byrow=TRUE))
par(xaxt="n",yaxt="n",mar=c(1.5,1,2,1))
for(t in 1:n.years){
  year.tmp=(fullD$Year==years[t])
  cover.tmp=fullD[year.tmp,'Cover']
  cover.mean.vec[t]=mean(cover.tmp,na.rm=TRUE)
  cover.center.tmp=cover.tmp-cover.mean.vec[t]
  cover.mean.center.vec[t]=mean(cover.center.tmp,na.rm=TRUE)
  lon.locs=sort(unique(fullD[year.tmp,'Lon']))
  lat.locs=sort(unique(fullD[year.tmp,'Lat']))
  n.lon=length(lon.locs)
  n.lat=length(lat.locs)
  image(matrix(cover.center.tmp,n.lon,n.lat),x=lon.locs,y=lat.locs,col=gray.colors.rev(100),asp=TRUE,main=years[t])
}
dev.off()

####
####  Plot Ratio Maps 
####

ratio.mat=matrix(0,n.sites,n.years)
png(file="sage_ratio_maps.png",width=10,height=5,units="in",res=200)
layout(matrix(1:n.years,4,7,byrow=TRUE))
par(xaxt="n",yaxt="n",mar=c(1.5,1,2,1))
for(t in 2:n.years){
  year.1.tmp=(fullD$Year==years[t-1])
  year.2.tmp=(fullD$Year==years[t])
  cover.1.tmp=fullD[year.1.tmp,'Cover']
  cover.2.tmp=fullD[year.2.tmp,'Cover']
  cover.ratio.tmp=cover.2.tmp/cover.1.tmp
  cover.ratio.tmp[cover.ratio.tmp>1000000]=NA
  ratio.mat[,t]=cover.ratio.tmp
  lon.locs=sort(unique(fullD[year.tmp,'Lon']))
  lat.locs=sort(unique(fullD[year.tmp,'Lat']))
  n.lon=length(lon.locs)
  n.lat=length(lat.locs)
  image(matrix(cover.ratio.tmp,n.lon,n.lat),x=lon.locs,y=lat.locs,col=gray.colors.rev(100),asp=TRUE,main=c(years[t-1],years[t]))
}
dev.off()

plot(years[-1],apply(ratio.mat[,-1],2,mean,na.rm=TRUE),type="o")
abline(h=1,col=8)
