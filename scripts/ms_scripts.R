###########################################################################################
##  ms_scripts.R: makes several figures for the manuscript based on the fitted sage      ##
##  abundance model. In particular, this script includes the simulations for Figure 4,   ##
##  which are simulations using equilibrium climate at 2098 and mean parameter estimates ##
##  The script also generates Figs. A1, F1, and 3.                                       ##
###########################################################################################

# Clear workspace 
rm(list=ls(all.names = TRUE))

####
####  Set some file paths, etc. ------------------------------------------------
####
datapath <- "../data/"
knotpath <- "../results/"



####
####  Load required libraries --------------------------------------------------
####
library(rstan)
library(reshape2)
library(ggplot2)
library(plyr)
library(sageAbundance)
library(gridExtra)
library(RColorBrewer)



####
####  Get data -----------------------------------------------------------------
####
fullD <- read.csv(paste0(datapath,"wy_sagecover_subset_noNA.csv"))
if(length(which(is.na(fullD$Cover))) > 0) stop("data contains NA values")

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

climD <- read.csv(paste0(datapath,
                        "/climate/DAYMET/FormattedClimate_WY_SA1.csv"))



####
####  Load knot data and MCMC results ------------------------------------------
####
load(paste0(knotpath,"Knot_cell_distances_subset.Rdata"))
K <- K.data$K
outs <- readRDS("../results/poissonSage_randYear_mcmc.RDS")



####
####  Fit colonization logistic ------------------------------------------------
####
colD <- growD[which(growD$CoverLag == 0), ]
colD$colonizes <- ifelse(colD$Cover==0, 0, 1)
col.fit <- glm(colonizes ~ 1, data=colD, family = "binomial")
col.intercept <- as.numeric(coef(col.fit))
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
avg.new.cover <- round(mean(colD[which(colD$Cover>0),"Cover"]),0)



####
####  Set up plotting themes ---------------------------------------------------
####
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=16),text=element_text(size=16),
                legend.text=element_text(size=16))
tmp.theme2=theme(strip.text=element_text(face="bold"),
                axis.title=element_text(size=16),text=element_text(size=16),
                legend.text=element_text(size=16))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))



####
####  Plot all years in space --------------------------------------------------
####
png("../results/all_years_percCover.png", width = 8.5, height = 5, 
    units="in", res = 100)
ggplot(growD, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=Cover, fill=Cover))+
  facet_wrap("Year")+
  scale_fill_gradientn(colours=myPalette(100))+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))
dev.off()



####
####  Plot spatial field and avg cover together --------------------------------
####
alpha <- outs[grep("alpha", outs$Parameter),]
alpha_means <- ddply(alpha, .(Parameter), summarise,
                     mean = mean(value))
eta <- K%*%alpha_means$mean
oneyear <- subset(growD, Year==1985)
oneyear$eta <- eta

png("../results/spatialfield.png", width = 8.5, height = 8.5, 
    units="in", res = 100)
etaPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
spfield_plot <- ggplot(oneyear, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=eta, fill=eta))+
  scale_fill_gradientn(limits=c(-0.55, 0.55), colours=etaPalette(100), 
                       name=expression(eta))+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle(bquote("Mean posterior spatial field ("*eta*")"))

avg_cover <- ddply(growD, .(Lon, Lat), summarise,
                   Cover = mean(Cover))
avg_cover$scaled_cov <- scale(avg_cover$Cover, center = TRUE, scale = TRUE)
avgcover_plot <- ggplot(avg_cover, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=scaled_cov, fill=scaled_cov))+
  scale_fill_gradientn(limits=c(-3, 3), colours=etaPalette(100), 
                       name="Cover \ndeviation")+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("Deviation from mean percent cover")

eta_avgcover_plot <- grid.arrange(avgcover_plot, spfield_plot, nrow=2)
dev.off()


####
####  Compare spatial field to landscape features ------------------------------
####
# slope <- raster("/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/ancillary_data/ls3731_slope_subset.img")
# slopeDf <- as.data.frame(rasterToPoints(slope))
# eta_slope <- merge(oneyear, slopeDf, by.x=c("Lon", "Lat"), by.y=c("x", "y"))
# ggplot(eta_slope, aes(x=ls3731_slope_subset, y=eta))+
#   geom_point()+
#   stat_smooth(se=FALSE, color="black", method="gam", formula=y ~ s(x, bs = "cs"))+
#   theme_bw()


####
####  Plot climate covariate posterior densities -------------------------------
####
clim_vars <- c("beta.1.", "beta.2.", "beta.3.", "beta.4.", "beta.5.")
clim_names <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
clim_posts <- subset(outs, Parameter %in% clim_vars)
quants <- ddply(clim_posts, .(Parameter), summarise,
              lower = quantile(value, .05),
              upper = quantile(value, .95),
              average = mean(value))
quants$clim_names <- clim_names
quants <- quants[with(quants, order(-average)), ]
quants$rank <- c(1:nrow(quants))
quants$col <- NA

for(i in 1:nrow(quants)){
  if(sign(quants$lower[i]) == sign(quants$upper[i]) & quants$average[i] > 0 ){
    quants$col[i] <- "pos"
  }
  if(sign(quants$lower[i]) == sign(quants$upper[i]) & quants$average[i] < 0 ){
    quants$col[i] <- "neg"
  }
  
  if(sign(quants$lower[i]) != sign(quants$upper[i]) ){
    quants$col[i] <- "none"
  }
}

clim_posts <- merge(clim_posts, quants, by="Parameter")

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#CC79A7", "#0072B2", "#D55E00", "#CC79A7")
my_labeller <- function(variable, value){
  return(as.character(quants$clim_names))
} 
library(ggthemes)
png("../results/post_climate_covariates.png", height=4, width=4, units="in", res=100)
ggplot(clim_posts, aes(x=value, fill=col))+
  geom_vline(aes(xintercept=0), linetype=2)+
  geom_density(alpha=0.5,size=0.5, col=NA, adjust=4)+
  scale_fill_manual(values=c("grey60", "steelblue", "red"))+
  ylab("Posterior density")+
  xlab("Standardized coefficient value")+
  facet_grid(rank~., labeller=my_labeller)+
  scale_y_continuous(breaks=NULL)+
  theme_few()+
  tmp.theme2+
  theme(strip.background=element_rect(fill="white"))+
  theme(strip.text.y = element_text(size = 8, angle = 270))+
  guides(fill=FALSE)
dev.off()



####
####  Run equilibrium population simulation ------------------------------------
####
time.steps <- 2000
burn.in <- 100
mean_params <- ddply(outs, .(Parameter), summarise,
                     value = mean(value))
yrint_ids <- grep("int_yr", mean_params$Parameter)
int_yrs <- mean_params[yrint_ids,"value"]
alphas <- mean_params[grep("alpha", mean_params$Parameter),"value"]
betas <- mean_params[grep("beta", mean_params$Parameter),"value"][2:6]
eta <- K%*%alphas
pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
clim_sim <- climD[climD$year %in% unique(growD$Year),]
X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 2:time.steps){
  int_now <- int_yrs[sample(c(1:length(int_yrs)), 1)]
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  dens.dep <- mean_params[mean_params$Parameter=="beta_mu","value"]*log(ex.mat[t-1,])
  tmp.mu <- int_now + dens.dep + sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  
  #Colonization
  zeros <- which(ex.mat[t-1,]==0)
  colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
  colonizer.cover <- colonizers*avg.new.cover
  tmp.out[zeros] <- colonizer.cover
  
  ex.mat[t,] <- tmp.out
}

# matplot(c(1:time.steps), ex.mat, type="l", col="grey")


####
####  Plot spatial equilibrium -------------------------------------------------
####
obs.equil <- ddply(growD, .(Lon, Lat), summarise,
                   cover = mean(Cover))
obs.equil <- obs.equil[with(obs.equil, order(-Lat, Lon)), ]
eq.pixels <- colMeans(ex.mat[(burn.in+1):time.steps,])
sp.equil <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
                       Lat=subset(growD, Year==1985)$Lat,
                       pred.cover=eq.pixels,
                       obs.cover=obs.equil$cover)
colnames(sp.equil)[3:4] <- c("Predicted", "Observed") 
sp.equil2 <- melt(sp.equil, id.vars = c("Lon", "Lat"))

##  Calculate prediction bias and precision
wide_equil <- dcast(sp.equil2, Lon+Lat~variable, value.var = "value")
wide_equil$bias <- with(wide_equil, Predicted-Observed)

ex2 <- ex.mat[(burn.in+1):time.steps,]
pred.var <- apply(ex2, 2, sd)
sp.predvar <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
                       Lat=subset(growD, Year==1985)$Lat,
                       pred.prec=1/pred.var)

tmp.theme3=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=12),text=element_text(size=12),
                legend.text=element_text(size=10))

png("../results/obs_pred_spatial_subset.png", width = 8, height=8, units = "in", res=300)
g1 <- ggplot(subset(sp.equil2, variable=="Observed"), aes(x=Lon, y=Lat))+
  geom_raster(aes(z=value, fill=value))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover", limit=c(0,25))+
  coord_equal()+
  tmp.theme3+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("A) Observed cover")

g2 <- ggplot(subset(sp.equil2, variable=="Predicted"), aes(x=Lon, y=Lat))+
  geom_raster(aes(z=value, fill=value))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover", limit=c(0,25))+
  coord_equal()+
  tmp.theme3+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("B) Predicted cover")


g3 <- ggplot(wide_equil, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=bias, fill=bias))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover", limits=c(-20, 20))+
  coord_equal()+
  tmp.theme3+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("C) Prediction bias")

g4 <- ggplot(sp.predvar, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=pred.prec, fill=pred.prec))+
  scale_fill_gradientn(colours=myPalette(200), name=expression("(% Cover)"^-2))+
  coord_equal()+
  tmp.theme3+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("D) Prediction precision")

gout <- grid.arrange(g1, g2, g3, g4, ncol=1, nrow=4)
dev.off()


####
####  Predict each year --------------------------------------------------------
####
mean_params <- ddply(outs, .(Parameter), summarise,
                     mean(value))
yrint_ids <- grep("int_yr", mean_params$Parameter)
int_yrs <- mean_params[yrint_ids,"value"]
alphas <- mean_params[grep("alpha", mean_params$Parameter),"..1"]
betas <- mean_params[grep("beta", mean_params$Parameter),"..1"][2:6]
eta <- K%*%alphas
time.steps <- length(unique(growD$Year))
yearid <- unique(growD$Year)
pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
clim_sim <- climD[climD$year %in% unique(growD$Year),]
X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
X_sim = scale(X_sim, center = TRUE, scale = TRUE)
for(t in 1:time.steps){
  int_now <- int_yrs[t]
  Xtmp <- X_sim[t,]
  lagcover <- growD[which(growD$Year==yearid[t]),"CoverLag"]
  dens.dep <- mean_params[mean_params$Parameter=="beta_mu","..1"]*log(lagcover)
  dens.dep[which(dens.dep==-Inf)] <- 0
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + dens.dep + sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  
  #Colonization
  zeros <- which(lagcover==0)
  colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
  colonizer.cover <- colonizers*avg.new.cover
  tmp.out[zeros] <- colonizer.cover
  
  ex.mat[t,] <- tmp.out
  
}

obs.cover <- growD[,c("Year","Cover","ID")]

ex.df <- as.data.frame(ex.mat)
colnames(ex.df) <- unique(obs.cover$ID)
ex.df$year <- yearid
ex.melt <- melt(ex.df, id.vars = "year")
ex.melt <- ex.melt[with(ex.melt, order(year)), ]
colnames(ex.melt) <- c("Year", "ID", "cover.pred")
error.df <- merge(obs.cover, ex.melt)
error.df$error <- with(error.df, Cover - cover.pred)
insamp.rmse <- sqrt(mean(error.df$error^2))
insamp.cor <- cor(error.df$Cover, error.df$cover.pred)
out <- data.frame(insamp_rmse = insamp.rmse,
                  insamp_cor = insamp.cor)
write.csv(out, "../results/insample_error_cor.csv")

