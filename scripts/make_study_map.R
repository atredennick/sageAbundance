##  Script to make study area map

rm(list=ls()) #wipe workspace clean

require(ggplot2)
require(ggthemes)
require(raster)
require(maptools)
require(shapefiles)
require(sp)
require(rgdal)
library(colorRamps)

#Source scale_bar function
source("/Users/atredenn/Desktop/ggplot_scale_bar.R")

#Get multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

####
####  Read in and format map data ----------------------------------------------
####
gpclibPermit()
US.shp <- readShapeSpatial("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/GIS/US_Shapefiles/Western_USt.shp")
US.shp@data$id <- rownames(US.shp@data)
US.pts <- fortify(US.shp, region="id")
WY.shp <- readShapeSpatial("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/GIS/US_Shapefiles/Wyoming.shp")
WY.shp@data$id <- rownames(WY.shp@data)
WY.pts <- fortify(WY.shp, region="id")



####
####  Read in a format sagebrush RS data ---------------------------------------
####
# 1.1-Bookkeeping to keep things in chronological order
setwd("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/sage/")
fileList <- as.matrix(list.files())
years <- c(seq(2000,2011,1), seq(1984,1999,1))
imgYears <- cbind(fileList, years) 
imgYears <- imgYears[order(imgYears[,2]),]

#Just get 1 year of data for plotting
shrub.ras <- raster(imgYears[1,1])
shrub.ras <- projectRaster(shrub.ras, crs=CRS("+proj=longlat"))
shrub.ras[shrub.ras>40] <- NA
shrub.ras[shrub.ras<0] <- NA
shrub.df <- data.frame(rasterToPoints(shrub.ras))
shrub.df$Cover <- round(shrub.df[,3])
shrub.center <- median(seq(1,dim(shrub.df)[1]))


setwd("/Users/atredenn/Repos/sageAbundance/scripts/")
sage_data <- read.csv("../data/wy_sagecover_subset_noNA.csv")
sage_oneyear <- subset(sage_data, Year==1985)
rm(sage_data)

ggplot()+
  geom_polygon(data=US.pts, aes(long,lat, group=group), color="grey50", 
               fill="white", lineend = "round")+
  geom_point(data=shrub.df[shrub.center,], aes(y=y, x=x), 
             color="black", fill="black", size=2, shape=22) +
  coord_equal() +
  guides(color=FALSE)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
  )+
  scaleBar(lon = -130, lat = 26, distanceLon = 500, distanceLat = 100, 
           distanceLegend = 200, dist.unit = "km")
ggsave("../docs/components/figure/us_wyo_studyarea_mark.png")

ggplot()+
  geom_raster(data=shrub.df, aes(y=y, x=x, fill=Cover))+
  theme_bw() +
  coord_equal() +
  scale_fill_gradient(high="white", low="black") +
  xlab("Longitude") + ylab("Latitude")+ 
  theme(
        panel.background=element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks=element_blank()
  )
ggsave("../docs/components/figure/oneyear_cover_example.png")





