## Devel scripts for manscript

# Clear workspace 
rm(list=ls(all=TRUE))

####
####  Set some file paths, etc.
####
datapath <- "/Users/atredenn/Dropbox/sageAbundance_data/"
knotpath <- "../results/"

####
####  Load required libraries
####
library(rstan)
library(reshape2)
library(ggplot2)
library(plyr)
library(sageAbundance)
library(gridExtra)
library(RColorBrewer)



####
####  Get data
####
fullD <- read.csv(paste0(datapath,"wy_sagecover_subset_noNA.csv"))
if(length(which(is.na(fullD$Cover))) > 0) stop("data contains NA values")

# Get data structure right
growD <- subset(fullD, Year>1984) # get rid of NA lagcover years
growD$Cover <- round(growD$Cover,0) # round for count-like data
growD$CoverLag <- round(growD$CoverLag,0) # round for count-like data

climD <- read.csv(paste(datapath,
                        "/studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv",
                        sep=""))

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=16),text=element_text(size=16),
                legend.text=element_text(size=16))
tmp.theme2=theme(strip.text=element_text(face="bold"),
                axis.title=element_text(size=16),text=element_text(size=16),
                legend.text=element_text(size=16))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#all years cover plot

png("../results/all_years_percCover.png", width = 8.5, height = 5, units="in", res = 100)
ggplot(growD, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=Cover, fill=Cover))+
  facet_wrap("Year")+
  scale_fill_gradientn(colours=myPalette(100))+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))
dev.off()

# Load knot data
load(paste0(knotpath,"Knot_cell_distances_subset.Rdata"))
K <- K.data$K

outs <- readRDS("../results/poissonSage_mcmc.RDS")

## Look at spatial field
alpha <- outs[grep("alpha", outs$Parameter),]
alpha_means <- ddply(alpha, .(Parameter), summarise,
                     mean = mean(value))
eta <- K%*%alpha_means$mean
oneyear <- subset(growD, Year==1985)
oneyear$eta <- eta



####
####  Plot spatial field and avg cover together
####
png("../results/spatialfield.png", width = 8.5, height = 8.5, units="in", res = 100)
etaPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
spfield_plot <- ggplot(oneyear, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=eta, fill=eta))+
  scale_fill_gradientn(limits=c(-0.55, 0.55), colours=etaPalette(100), name=expression(eta))+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle(bquote("Mean posterior spatial field ("*eta*")"))

avg_cover <- ddply(growD, .(Lon, Lat), summarise,
                   Cover = mean(Cover))
avg_cover$scaled_cov <- scale(avg_cover$Cover, center = TRUE, scale = TRUE)
avgcover_plot <- ggplot(avg_cover, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=scaled_cov, fill=scaled_cov))+
  scale_fill_gradientn(limits=c(-3, 3), colours=etaPalette(100), name="Cover \ndeviation")+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("Deviation from mean percent cover")

eta_avgcover_plot <- grid.arrange(avgcover_plot, spfield_plot, nrow=2)
dev.off()


####
####  Plot climate covariate posterior densities
####
clim_vars <- c("beta.1.", "beta.2.", "beta.3.", "beta.4.", "beta.5.")
clim_names <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
clim_posts <- subset(outs, Parameter %in% clim_vars)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#CC79A7", "#0072B2", "#D55E00", "#CC79A7")
png("../results/post_climate_covariates.png", height=5, width=8, units="in", res=100)
ggplot(clim_posts, aes(x=value, color=Parameter))+
  geom_vline(aes(xintercept=0), linetype=2)+
  geom_line(stat="density", size=1.5)+
  scale_color_manual(labels=clim_names, values=cbPalette, name="Covariate")+
  ylab("Posterior density")+
  xlab("Standardized coefficient value")+
  tmp.theme2+
  theme(strip.background=element_rect(fill="white"))
dev.off()



####
####  Run population simulation
####
time.steps <- 600
burn.in <- 100
mean_params <- ddply(outs, .(Parameter), summarise,
                     value = mean(value))
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
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","value"] + mean_params[mean_params$Parameter=="beta_mu","value"]*ex.mat[t-1,] + 
    sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
#   print(max(tmp.mu))
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
  ex.mat[t,] <- tmp.out
}

# matplot(c(1:time.steps), ex.mat, type="l", col="grey")
# # hist(y, freq = FALSE)
# plot(density(ex.mat, adjust = 4), col="red", lwd=2)
# abline(v = mean(growD$Cover))
# mean(ex.mat)


####
####  Plot spatial equilibrium
####
obs.equil <- ddply(growD, .(Lon, Lat), summarise,
                   cover = mean(Cover))
obs.equil <- obs.equil[with(obs.equil, order(-Lat, Lon)), ]
# obs.equil[which(obs.equil[,"cover"]>10), "cover"] <- NA
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

pred.var <- apply(ex.mat, 2, var)
sp.predvar <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
                       Lat=subset(growD, Year==1985)$Lat,
                       pred.prec=1/pred.var)

tmp.theme3=theme(axis.ticks = element_blank(), axis.text = element_blank(),
                strip.text=element_text(face="bold"),
                axis.title=element_text(size=12),text=element_text(size=12),
                legend.text=element_text(size=10))

png("../results/obs_pred_spatial_subset.png", width = 8, height=8, units = "in", res=100)
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
  scale_fill_gradientn(colours=myPalette(200), name="% Cover")+
  coord_equal()+
  tmp.theme3+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("C) Prediction bias")

g4 <- ggplot(sp.predvar, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=pred.prec, fill=pred.prec))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover")+
  coord_equal()+
  tmp.theme3+
  theme(strip.background=element_rect(fill="white"))+
  ggtitle("D) Prediction precision")

# g2 <- ggplot(sp.equil)+
#   geom_histogram(aes(x=Observed, y = ..density..), col="white", fill="grey45", binwidth=1)+
#   geom_line(stat="density", aes(x=Predicted, y= ..density..), 
#             col="black", size=1.5, linetype=2)+
#   xlab("Percent Cover")+
#   ylab("Estimated Kernel Density")+
#   theme_bw()
gout <- grid.arrange(g1, g2, g3, g4, ncol=1, nrow=4)
dev.off()


####
####  Predict each year
####
mean_params <- ddply(outs, .(Parameter), summarise,
                     mean(value))
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
  Xtmp <- X_sim[t,]
  lagcover <- growD[which(growD$Year==yearid[t]),"CoverLag"]
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*lagcover + 
    sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
  ex.mat[t,] <- tmp.out
}

obs.cover <- growD[,c("Year","Cover","ID")]
pred.cover <- melt(ex.mat)
plot(obs.cover$Cover, pred.cover$value, xlim=c(0,20), ylim=c(0,20))
abline(0,1, col="red")

brks <- 30
y <- hist(growD$Cover, freq=FALSE, breaks=brks)
yhat <- hist(pred.cover$value, freq=FALSE, breaks=brks)
hist_ydata <- data.frame(x = y$breaks, y=c(y$density,0))
hist_yhatdata <- data.frame(x = yhat$breaks, y=c(yhat$density,0))

ggplot()+
  geom_step(data=hist_ydata, aes(x=x, y=y), size=1.5)+
  geom_step(data=hist_yhatdata, aes(x=x, y=y), col="darkorange", size=1)+
  xlab("Cover (%)")+
  ylab("Density")
ggsave("../results/obs_pred_hist.png", width = 4, height = 4, units="in")





####
####  Run climate simulations
####
projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
                  "deltaPpt"=c(1.0894,1.0864,1.110),
                  "deltaTspr"=c(2.975,3.134,4.786))
mean_params <- ddply(outs, .(Parameter), summarise,
                     mean(value))
alphas <- mean_params[grep("alpha", mean_params$Parameter),"..1"]
betas <- mean_params[grep("beta", mean_params$Parameter),"..1"][2:6]
eta <- K%*%alphas
time.steps <- 600
burn.in <- 100

p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
clim_avg <- apply(X = p.climD, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = p.climD, MARGIN = 2, FUN = sd)



pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[1,2],dim(climD)[1],2)
p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[1,3],dim(climD)[1],2)
X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]

# Now scale based on perturbed or regular data, depending on scenario
X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
                  sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  #   print(max(tmp.mu))
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
  ex.mat[t,] <- tmp.out
}
rcp45 <- ex.mat[(burn.in+1):time.steps, ]


pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[2,2],dim(climD)[1],2)
p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[2,3],dim(climD)[1],2)
X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# Now scale based on perturbed or regular data, depending on scenario
X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
    sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  #   print(max(tmp.mu))
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
  ex.mat[t,] <- tmp.out
}
rcp60 <- ex.mat[(burn.in+1):time.steps, ]


pixels <- nrow(subset(growD, Year==1985))
ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
ex.mat[1,] <- 1
p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[3,2],dim(climD)[1],2)
p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[3,3],dim(climD)[1],2)
X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# Now scale based on perturbed or regular data, depending on scenario
X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
for(t in 2:time.steps){
  Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
  tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
    sum(betas*Xtmp)
  tmp.mu <- exp(tmp.mu + eta)
  #   print(max(tmp.mu))
  tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
  tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
  ex.mat[t,] <- tmp.out
}
rcp85 <- ex.mat[(burn.in+1):time.steps, ]


####
####  Plot climate change results
####
all.sims <- data.frame(pixel=c(1:length(eq.pixels)),
                       Baseline=eq.pixels,
                       RCP45=colMeans(rcp45),
                       RCP60=colMeans(rcp60),
                       RCP85=colMeans(rcp85))
all4plot <- melt(all.sims, id.vars = "pixel")




####
####  Plot spatial projections for RCP scenarios
####
proj.equil <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
                         Lat=subset(growD, Year==1985)$Lat,
                         CURRENT=eq.pixels,
                         RCP45=colMeans(rcp45),
                         RCP60=colMeans(rcp60),
                         RCP85=colMeans(rcp85))
proj.equil2 <- melt(proj.equil, id.vars = c("Lon", "Lat"))

# png("../results/climchange_small.png", width = 6, height=4.5, units = "in", res=150)
ggplot(proj.equil2, aes(x=Lon, y=Lat))+
  geom_raster(aes(z=value, fill=value))+
  scale_fill_gradientn(colours=myPalette(200), name="% Cover")+
  facet_wrap("variable", ncol=1)+
  coord_equal()+
  tmp.theme+
  theme(strip.background=element_rect(fill="white"))
ggsave("../results/clim_change_mean_spatial.png", height=8, width=4)

ggplot(all4plot)+
  geom_line(stat="density", aes(x=value, y=..density.., col=variable), size=1)+
  scale_color_manual(values=c("dodgerblue","tan","coral","darkred"), name="Scenario") +
  xlab("Equilibrium Percent Cover")+
  ylab("Estimated Kernel Density")
ggsave("../results/clim_change_densities.png", height=4, width=5, dpi = 200)
# gout <- grid.arrange(g1, g2, ncol=1, nrow=2)
# dev.off()





# 
# 
# ###################################################
# ###################################################
# ###################################################
# ####
# ####  Run population model with parameter uncertainty
# ####
# colnames(outs)[which(colnames(outs)=="chain")] <- "Chain"
# 
# time.steps <- 10
# pixels <- nrow(subset(growD, Year==1985))
# clim_sim <- climD[climD$year %in% unique(growD$Year),]
# X_sim = clim_sim[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# X_sim = scale(X_sim, center = TRUE, scale = TRUE)
# 
# alphasd <- outs[grep("alpha", outs$Parameter),]
# climeffs <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")
# nchains <- 3
# niters <- length(unique(outs$mcmc_iter))
# totsims <- nchains*niters
# ex.arr <- array(NA, dim=c(totsims, time.steps, pixels))
# ex.arr[,1,] <- 1
# pb <- txtProgressBar(min=1, max=totsims, char="+", style=3, width=65)
# counter <- 1
# for(i in 1:nchains){
#   for(j in 1:niters){
#     chain <- i
#     iter <- j
#     int_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="int_mu")[,"value"])
#     beta_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="beta_mu")[,"value"])
#     betas <- subset(outs, Chain==chain & mcmc_iter==iter & Parameter%in%climeffs)[,"value"]
#     alphas <- subset(alphasd, Chain==chain & mcmc_iter==iter)[,"value"]
#     eta <- K%*%as.numeric(unlist(alphas))
#     for(t in 2:time.steps){
#       Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#       tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
#         sum(betas*Xtmp)
#       tmp.mu <- exp(tmp.mu + eta)
#       #   print(max(tmp.mu))
#       tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#       tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
#       ex.arr[counter,t,] <- tmp.out
#     }
#     counter <- counter+1
#     setTxtProgressBar(pb, counter)
#   }
# }
# 
# y <- hist(growD$Cover, freq=FALSE)
# yhat <- hist(ex.arr, freq=FALSE)
# hist_ydata <- data.frame(x = y$breaks, y=c(y$density,0))
# hist_yhatdata <- data.frame(x = yhat$breaks, y=c(yhat$density,0))
# # pdf("~/Desktop/test_hist.pdf", width = 4, height = 4)
# ggplot()+
#   geom_step(data=hist_ydata, aes(x=x, y=y), size=1.5)+
#   geom_step(data=hist_yhatdata, aes(x=x, y=y), col="darkorange", size=1)+
#   xlab("Cover (%)")+
#   ylab("Density")+
#   theme_bw()
# # dev.off()
# 
# 
# ####
# ####  RCP 4.5; parameter vary
# ####
# projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
#                   "deltaPpt"=c(1.0894,1.0864,1.110),
#                   "deltaTspr"=c(2.975,3.134,4.786))
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# clim_avg <- apply(X = p.climD, MARGIN = 2, FUN = mean)
# clim_sd <- apply(X = p.climD, MARGIN = 2, FUN = sd)
# p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[1,2],dim(climD)[1],2)
# p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[1,3],dim(climD)[1],2)
# X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# 
# # Now scale based on perturbed or regular data, depending on scenario
# X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
# X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
# X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
# X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
# X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
# 
# alphasd <- outs[grep("alpha", outs$Parameter),]
# climeffs <- c("beta.1.", "beta.2.", "beta.3.", "beta.4.", "beta.5.")
# nchains <- 3
# niters <- length(unique(outs$mcmc_iter))
# totsims <- nchains*niters
# pixels <- nrow(subset(growD, Year==1985))
# ex.arr <- array(NA, dim=c(totsims, time.steps, pixels))
# ex.arr[,1,] <- 1
# pb <- txtProgressBar(min=1, max=totsims, char="+", style=3, width=65)
# counter <- 1
# for(i in 1:nchains){
#   for(j in 1:niters){
#     chain <- i
#     iter <- j
#     int_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="int_mu")[,"value"])
#     beta_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="beta_mu")[,"value"])
#     betas <- subset(outs, Chain==chain & mcmc_iter==iter & Parameter%in%climeffs)[,"value"]
#     betas <- as.numeric(unlist(betas))
#     alphas <- subset(alphasd, Chain==chain & mcmc_iter==iter)[,"value"]
#     eta <- K%*%as.numeric(unlist(alphas))
#     for(t in 2:time.steps){
#       Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#       tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
#         sum(betas*Xtmp)
#       tmp.mu <- exp(tmp.mu + eta)
#       #   print(max(tmp.mu))
#       tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#       tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
#       ex.arr[counter,t,] <- tmp.out
#     }
#     counter <- counter+1
#     setTxtProgressBar(pb, counter)
#   }
# }
# rcp45 <- ex.arr
# rcp45a <- alply(rcp45,1)
# 
# 
# ####
# ####  RCP 6.0; parameter vary
# ####
# projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
#                   "deltaPpt"=c(1.0894,1.0864,1.110),
#                   "deltaTspr"=c(2.975,3.134,4.786))
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# clim_avg <- apply(X = p.climD, MARGIN = 2, FUN = mean)
# clim_sd <- apply(X = p.climD, MARGIN = 2, FUN = sd)
# p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[2,2],dim(climD)[1],2)
# p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[2,3],dim(climD)[1],2)
# X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# 
# # Now scale based on perturbed or regular data, depending on scenario
# X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
# X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
# X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
# X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
# X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
# 
# alphasd <- outs[grep("alpha", outs$Parameter),]
# nchains <- 3
# niters <- length(unique(outs$mcmc_iter))
# totsims <- nchains*niters
# pixels <- nrow(subset(growD, Year==1985))
# ex.arr <- array(NA, dim=c(totsims, time.steps, pixels))
# ex.arr[,1,] <- 1
# pb <- txtProgressBar(min=1, max=totsims, char="+", style=3, width=65)
# counter <- 1
# for(i in 1:nchains){
#   for(j in 1:niters){
#     chain <- i
#     iter <- j
#     int_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="int_mu")[,"value"])
#     beta_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="beta_mu")[,"value"])
#     betas <- subset(outs, Chain==chain & mcmc_iter==iter & Parameter%in%climeffs)[,"value"]
#     betas <- as.numeric(unlist(betas))
#     alphas <- subset(alphasd, Chain==chain & mcmc_iter==iter)[,"value"]
#     eta <- K%*%as.numeric(unlist(alphas))
#     for(t in 2:time.steps){
#       Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#       tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
#         sum(betas*Xtmp)
#       tmp.mu <- exp(tmp.mu + eta)
#       #   print(max(tmp.mu))
#       tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#       tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
#       ex.arr[counter,t,] <- tmp.out
#     }
#     counter <- counter+1
#     setTxtProgressBar(pb, counter)
#   }
# }
# rcp60 <- ex.arr
# rcp60a <- alply(rcp60,1)
# 
# ####
# ####  RCP 8.5; parameter vary
# ####
# projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
#                   "deltaPpt"=c(1.0894,1.0864,1.110),
#                   "deltaTspr"=c(2.975,3.134,4.786))
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# clim_avg <- apply(X = p.climD, MARGIN = 2, FUN = mean)
# clim_sd <- apply(X = p.climD, MARGIN = 2, FUN = sd)
# p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[3,2],dim(climD)[1],2)
# p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[3,3],dim(climD)[1],2)
# X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# 
# # Now scale based on perturbed or regular data, depending on scenario
# X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
# X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
# X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
# X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
# X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
# 
# alphasd <- outs[grep("alpha", outs$Parameter),]
# nchains <- 3
# niters <- length(unique(outs$mcmc_iter))
# totsims <- nchains*niters
# pixels <- nrow(subset(growD, Year==1985))
# ex.arr <- array(NA, dim=c(totsims, time.steps, pixels))
# ex.arr[,1,] <- 1
# pb <- txtProgressBar(min=1, max=totsims, char="+", style=3, width=65)
# counter <- 1
# for(i in 1:nchains){
#   for(j in 1:niters){
#     chain <- i
#     iter <- j
#     int_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="int_mu")[,"value"])
#     beta_mu <- as.numeric(subset(outs, Chain==chain & mcmc_iter==iter & Parameter=="beta_mu")[,"value"])
#     betas <- subset(outs, Chain==chain & mcmc_iter==iter & Parameter%in%climeffs)[,"value"]
#     betas <- as.numeric(unlist(betas))
#     alphas <- subset(alphasd, Chain==chain & mcmc_iter==iter)[,"value"]
#     eta <- K%*%as.numeric(unlist(alphas))
#     for(t in 2:time.steps){
#       Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#       tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + mean_params[mean_params$Parameter=="beta_mu","..1"]*ex.mat[t-1,] + 
#         sum(betas*Xtmp)
#       tmp.mu <- exp(tmp.mu + eta)
#       #   print(max(tmp.mu))
#       tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#       tmp.out[which(tmp.out>100)] <- mean(ex.mat[t-1,])
#       ex.arr[counter,t,] <- tmp.out
#     }
#     counter <- counter+1
#     setTxtProgressBar(pb, counter)
#   }
# }
# rcp85 <- ex.arr
# rcp85a <- alply(rcp85,1)
# 
# 
# ####
# ####  Spatial plot of forecast diffs; cover only for >0.6 prob of change
# ####
# calcdiffs <- function(ypred.list, yobs, limit=TRUE, limit.val=0.6){
#   if(dim(ypred.list[[1]])[2] != length(yobs)) stop("dimension mismatch between ypred and yobs")
#   n <- length(ypred.list)
#   tmpdiff <- matrix(ncol=dim(ypred.list[[1]])[2], nrow=n)
#   for(i in 1:n){
#     tmpdiff[i,] <- colMeans(ypred.list[[i]]) - yobs
#   }
#   
#   if(limit==TRUE){
#     out <- numeric(dim(tmpdiff)[2])
#     for(i in 1:length(out)){
#       tmp <- abs(ecdf(tmpdiff[,i])(0) - (1-ecdf(tmpdiff[,i])(0)))
#       if(tmp > limit.val) out[i] <- mean(tmpdiff[,i])
#       if(tmp <= limit.val) out[i] <- 0
#     }# end pixel loop
#   }# end limit if/then
#   return(out)
# }# end function
# 
# rcp45.diffs <- calcdiffs(rcp45a, yobs = obs.equil$cover, limit = TRUE, limit.val = 0.95)
# rcp60.diffs <- calcdiffs(rcp60a, yobs = obs.equil$cover, limit = TRUE, limit.val = 0.95)
# rcp85.diffs <- calcdiffs(rcp85a, yobs = obs.equil$cover, limit = TRUE, limit.val = 0.95)
# proj.diff <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
#                         Lat=subset(growD, Year==1985)$Lat,
#                         RCP45=rcp45.diffs,
#                         RCP60=rcp60.diffs,
#                         RCP85=rcp85.diffs)
# 
# proj.diff2 <- melt(proj.diff, id.vars = c("Lon", "Lat"))
# proj.diff2[which(proj.diff2[,"value"]< (-10)), "value"] <- NA
# max(proj.diff2$value, na.rm = TRUE)
# min(proj.diff2$value, na.rm = TRUE)
# myPalette2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
# ggplot(proj.diff2, aes(x=Lon, y=Lat))+
#   geom_raster(aes(z=value, fill=value))+
#   scale_fill_gradientn(colours=myPalette2(200), name="Cover\nDifference (%)", limit=c(-20,20))+
#   facet_wrap("variable", nrow=3)+
#   coord_equal()+
#   tmp.theme+
#   theme(strip.background=element_rect(fill="white"))
# 
# 
# ggplot(proj.diff2, aes(x=value, color=variable))+
#   geom_density()
