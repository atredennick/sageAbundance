
library(reshape2)
library(raster)
library(sp)
library(rgeos)
library(Matrix)

####
#### Bring in data
####
climD<-read.csv("../../studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv")
rawD<-read.csv("../../cover_structure/WY_SAGECoverData_V2small.csv")


# merge in climate data 
fullD <- merge(rawD,climD,by.x="Year", by.y="year",all.x=T)

shrub.ras <- raster("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/sage/ls3731_00sage_subset_est.img")
shrub.ras <- aggregate(shrub.ras, fact=10, fun="mean", na.rm=TRUE)

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
#### Generate knots and plot
####
Coords=coordinates(dat)
x.min=min(Coords[,1])
x.max=max(Coords[,1])
y.min=min(Coords[,2])
y.max=max(Coords[,2])

splits <- 2
X=x.min+(x.max-x.min)/splits*c(0:splits)
Y=y.min+(y.max-y.min)/splits*c(splits:0)
XY=expand.grid(x=X,y=Y)

Knots=SpatialPoints(coords=XY,proj4string=CRS(proj4string(shrub.ras)))
obsPts=SpatialPoints(coords=Coords, proj4string=CRS(proj4string(shrub.ras)))
Distances=gDistance(Knots,obsPts,byid=TRUE)
Distances=apply(Distances,2,'min')
my.buffer=15000
Which.include=which(Distances<my.buffer)
# save(Knots,file="BOSS_Knots_SP.Rda")
Knot.cell.distances=gDistance(Knots[Which.include,],obsPts,byid=TRUE)
diff.x=(x.max-x.min)/splits #normally 6
diff.y=(y.max-y.min)/splits #normally 6
# sigma=(diff.x+diff.y)/2
sigma=250

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
save(K.data,file="Knot_cell_distances_smallSet.Rdata")

# png('SAGE_Grid_wKnots_smallSet.png', width=4, height=4, units = "in", res=200)
image(resMat, x=lon.locs,y=lat.locs,asp=TRUE,col=gray.colors.rev(100), 
      ylab="Latitude", xlab="Longitude")
points(Knots[Which.include,],,col='red',pch=20)
# dev.off()

Knots=Knots[Which.include,]
# save(Knots,file="SAGE_Knots_SP_smallSet.Rda")




