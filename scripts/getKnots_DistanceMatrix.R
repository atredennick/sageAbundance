
library(reshape2)
library(raster)
library(sp)
library(rgeos)
library(Matrix)
library(gstat)
library(plyr)
library(ggplot2)


####
#### Bring in data
####
setwd("/Users/atredenn/Dropbox/sageAbundance_data/")
studyArea="1"

ifelse(studyArea=="1",
       climD<-read.csv("./studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv"),
       climD<-read.csv("./studyarea2/climate/DAYMET/FormattedClimate_WY_SA2.csv"))
ifelse(studyArea=="1",
       rawD<-read.csv("./WY_SAGECoverData_V2check.csv"),
       rawD<-read.csv("./WY_SAGECoverData_SA2.csv"))

# merge in climate data 
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)

shrub.ras <- raster("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/sage/ls3731_00sage_subset_est.img")

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
subD <- subset(fullD, Year==years[1])
subC <- dcast(subD[,c('Lon', 'Lat', 'Cover')], formula = Lon~Lat)
resMat <- as.matrix(subC[, 2:ncol(subC)])
lon.locs=sort(unique(subD$Lon))
lat.locs=sort(unique(subD$Lat))
dat <- subset(subD, select = c("Lon", "Lat", "Cover"))
coordinates(dat) <- c("Lon", "Lat")
proj4string(dat)
image(resMat, x=lon.locs,y=lat.locs,asp=TRUE,col=gray.colors.rev(100))


####
#### Run simple Poisson GLM to get residual surface
####
modelD <- subset(fullD, Year>1984)
modelD <- modelD[complete.cases(modelD),]
outp <- glm(Cover ~ CoverLag+ppt2+ppt1+pptLag+TmeanSpr1+TmeanSpr2,
            family = "poisson", data = modelD)
summary(outp)
modelD$resids <- resid(outp)


####
####  Construct variograms to estimate range parameter for knots
####
resD <- ddply(modelD, .(Lat, Lon), summarise,
                  resid = mean(resids))
lon.locs=sort(unique(resD$Lon))
lat.locs=sort(unique(resD$Lat))
dat <- subset(resD, select = c("Lon", "Lat", "resid"))
coordinates(dat) <- c("Lon", "Lat")
varMod <- variogram(resid~1, data=dat, width=50)

png("avgresidual_variogram_withrange.png", width = 5, height = 5, 
    units = "in", res=200)
plot(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), col="dodgerblue", 
     type="l", main="variogram of model residuals", xlab="distance (m)",
     ylab="semivariance")
points(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), col="dodgerblue", pch=21, bg="white")
abline(v = 500/3, lwd=3, lty=2)
dev.off()

range_parameter <- round(500/3)

# par(xaxt="s",yaxt="s",mar= c(2, 2, 4, 2), mfrow=c(4,7))
# years <- unique(modelD$Year)
# for(t in years){
#   resD <- subset(modelD, Year==t)
#   lon.locs=sort(unique(resD$Lon))
#   lat.locs=sort(unique(resD$Lat))
#   dat <- subset(resD, select = c("Lon", "Lat", "resids"))
# #   dat <- dat[which(is.na(dat$Cover)==FALSE),]
#   coordinates(dat) <- c("Lon", "Lat")
#   varMod <- variogram(resids~1, data=dat, width=50)
#   plot(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), col="dodgerblue", type="l", main=years[t])
#   points(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), col="dodgerblue", pch=21, bg="white")
# }


####
#### Generate knots and plot
####
dat <- subset(resD, select = c("Lon", "Lat", "resid"))
coordinates(dat) <- c("Lon", "Lat")
Coords=coordinates(dat)
x.min=min(Coords[,1])
x.max=max(Coords[,1])
y.min=min(Coords[,2])
y.max=max(Coords[,2])

splits <- 20
X=x.min+(x.max-x.min)/splits*c(0:splits)
Y=y.min+(y.max-y.min)/splits*c(splits:0)
XY=expand.grid(x=X,y=Y)

Knots=SpatialPoints(coords=XY,proj4string=CRS(proj4string(shrub.ras)))
obsPts=SpatialPoints(coords=Coords, proj4string=CRS(proj4string(shrub.ras)))
Distances=gDistance(Knots,obsPts,byid=TRUE)
Distances=apply(Distances,2,'min')
my.buffer=150000
Which.include=which(Distances<my.buffer)
# save(Knots,file="BOSS_Knots_SP.Rda")
Knot.cell.distances=gDistance(Knots[Which.include,],obsPts,byid=TRUE)
diff.x=(x.max-x.min)/splits #normally 6
diff.y=(y.max-y.min)/splits #normally 6
test=(diff.x+diff.y)/2
sigma=range_parameter

#plot knot distances
Knot.distances=gDistance(Knots[Which.include,],Knots[Which.include,],byid=TRUE)
m <- melt(Knot.distances)
ggplot(data=m, aes(x=Var1, y=Var2))+
  geom_raster(aes(z=value, fill=value))

setwd("~/Repos/sageAbundance/R")
source("Conn_util_funcs.R")
Knot.Adj=rect_adj(splits+1,splits+1)
Knot.Adj=Knot.Adj[Which.include,Which.include]
Q.knot=-Knot.Adj
diag(Q.knot)=apply(Knot.Adj,2,'sum')
Q.knot=Matrix(Q.knot)

K=((2*pi*sigma)^-1)*exp(-Knot.cell.distances/(2*sigma))
# K=dnorm(Knot.cell.distances,0,sigma)
K=K/apply(K,1,'sum')
K.data=list(K=K,Q.knot=Q.knot)
save(K.data,file="Knot_cell_distances.Rdata")

#plot Conn's K
# load("/Users/atredenn/Desktop/STabundance-0.91/STabundance/data/Knot_cell_distances.rda")
# m <- melt(K.data$K)
# ggplot(data=m, aes(x=Var2, y=Var1))+
#   geom_raster(aes(z=value, fill=value))



library(plyr)
avgD <- ddply(fullD, .(Lon, Lat), summarise,
              cover = mean(Cover))
knotsD <- as.data.frame(Knots)

library(RColorBrewer)
myPalette <- colorRampPalette(brewer.pal(9, "Greens"))
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=14),text=element_text(size=16),
                legend.text=element_text(size=12), legend.title=element_text(size=16))

setwd("/Users/atredenn/Dropbox/sageAbundance_data/")
png('SAGE_Grid_wKnots.png', width=6, height=6, units = "in", res=200)
g <- ggplot()+
  geom_raster(data=avgD, aes(x=Lon, y=Lat, z=cover, fill=cover))+
  geom_point(data=knotsD, aes(x,y), size=1.5, color="white")+
  geom_point(data=knotsD, aes(x,y), size=1, color="black")+
  scale_fill_gradientn(colours=myPalette(100), name="Percent \nCover")+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))+
  coord_equal()+
  xlab("Longitude")+
  ylab("Latitude")
print(g)
dev.off()

# png('SAGE_Grid_wKnots.png', width=4, height=4, units = "in", res=200)
# image(resMat, x=lon.locs,y=lat.locs,asp=TRUE,col=gray.colors.rev(100), ylab="Latitude", xlab="Longitude")
# points(Knots[Which.include,],,col='red',pch=20, add=TRUE)
# dev.off()

Knots=Knots[Which.include,]
save(Knots,file="SAGE_Knots_SP.Rda")




