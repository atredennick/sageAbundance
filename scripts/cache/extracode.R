# pred.cover <- melt(ex.mat)
# plot(obs.cover$Cover, pred.cover$value, xlim=c(0,20), ylim=c(0,20))
# abline(0,1, col="red")
# brks <- 30
# y <- hist(growD$Cover, freq=FALSE, breaks=brks)
# yhat <- hist(pred.cover$value, freq=FALSE, breaks=brks)
# hist_ydata <- data.frame(x = y$breaks, y=c(y$density,0))
# hist_yhatdata <- data.frame(x = yhat$breaks, y=c(yhat$density,0))
# ggplot()+
#   geom_step(data=hist_ydata, aes(x=x, y=y), size=1.5)+
#   geom_step(data=hist_yhatdata, aes(x=x, y=y), col="darkorange", size=1)+
#   xlab("Cover (%)")+
#   ylab("Density")
# ggsave("../results/obs_pred_hist.png", width = 4, height = 4, units="in")





####
####  Run climate simulations --------------------------------------------------
####
##  This may take a while...
##  Runs equilibrium simulations for each model and RCP scenario

# Read in climate projected changes
# ppt_projs <- subset(read.csv("../data/precipitation_projections_bymodel.csv"), 
#                     scenario!="rcp26" & season=="fall2spr")
# temp_projs <- subset(read.csv("../data/temperature_projections_bymodel.csv"),
#                      scenario!="rcp26" & season=="spring")
# projC<-data.frame("scenario"=c("rcp45","rcp60","rcp85"),
#                   "deltaPpt"=1+ppt_projs$change,
#                   "deltaTspr"=temp_projs$change,
#                   "model"=ppt_projs$model)
# mean_params <- ddply(outs, .(Parameter), summarise,
#                      mean(value))
# alphas <- mean_params[grep("alpha", mean_params$Parameter),"..1"]
# betas <- mean_params[grep("beta", mean_params$Parameter),"..1"][2:6]
# eta <- K%*%alphas
# time.steps <- 2000
# burn.in <- 100
# 
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# clim_avg <- apply(X = p.climD, MARGIN = 2, FUN = mean)
# clim_sd <- apply(X = p.climD, MARGIN = 2, FUN = sd)
# 
# 
# 
# pixels <- nrow(subset(growD, Year==1985))
# ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
# ex.mat[1,] <- 10
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
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
# for(t in 2:time.steps){
#   Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#   dens.dep <- mean_params[mean_params$Parameter=="beta_mu","..1"]*log(ex.mat[t-1,])
#   dens.dep[which(dens.dep==-Inf)] <- 0
#   tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + dens.dep + sum(betas*Xtmp)
#   tmp.mu <- exp(tmp.mu + eta)
#   tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#   
#   #Colonization
#   zeros <- which(ex.mat[t-1,]==0)
#   colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
#   colonizer.cover <- colonizers*avg.new.cover
#   tmp.out[zeros] <- colonizer.cover
#   
#   ex.mat[t,] <- tmp.out
# }
# rcp45 <- ex.mat[(burn.in+1):time.steps, ]
# 
# 
# pixels <- nrow(subset(growD, Year==1985))
# ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
# ex.mat[1,] <- 10
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[2,2],dim(climD)[1],2)
# p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[2,3],dim(climD)[1],2)
# X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# # Now scale based on perturbed or regular data, depending on scenario
# X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
# X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
# X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
# X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
# X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
# for(t in 2:time.steps){
#   Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#   dens.dep <- mean_params[mean_params$Parameter=="beta_mu","..1"]*log(ex.mat[t-1,])
#   dens.dep[which(dens.dep==-Inf)] <- 0
#   tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + dens.dep + sum(betas*Xtmp)
#   tmp.mu <- exp(tmp.mu + eta)
#   tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#   
#   #Colonization
#   zeros <- which(ex.mat[t-1,]==0)
#   colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
#   colonizer.cover <- colonizers*avg.new.cover
#   tmp.out[zeros] <- colonizer.cover
#   
#   ex.mat[t,] <- tmp.out
# }
# rcp60 <- ex.mat[(burn.in+1):time.steps, ]
# 
# 
# pixels <- nrow(subset(growD, Year==1985))
# ex.mat <- matrix(NA,nrow=time.steps,ncol=pixels)
# ex.mat[1,] <- 10
# p.climD<-climD[climD$year %in% unique(growD$Year),c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# p.climD[,c(2:3)]<-p.climD[,c(2:3)]*matrix(projC[3,2],dim(climD)[1],2)
# p.climD[,c(4:5)]<-p.climD[,c(4:5)]+matrix(projC[3,3],dim(climD)[1],2)
# X_sim = p.climD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
# # Now scale based on perturbed or regular data, depending on scenario
# X_sim["pptLag"] <- (X_sim["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
# X_sim["ppt1"] <- (X_sim["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
# X_sim["ppt2"] <- (X_sim["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
# X_sim["TmeanSpr1"] <- (X_sim["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
# X_sim["TmeanSpr2"] <- (X_sim["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]
# for(t in 2:time.steps){
#   Xtmp <- X_sim[sample(c(1:nrow(X_sim)), 1),]
#   dens.dep <- mean_params[mean_params$Parameter=="beta_mu","..1"]*log(ex.mat[t-1,])
#   dens.dep[which(dens.dep==-Inf)] <- 0
#   tmp.mu <- mean_params[mean_params$Parameter=="int_mu","..1"] + dens.dep + sum(betas*Xtmp)
#   tmp.mu <- exp(tmp.mu + eta)
#   tmp.out <- rpois(ncol(ex.mat), lambda = tmp.mu)
#   
#   #Colonization
#   zeros <- which(ex.mat[t-1,]==0)
#   colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
#   colonizer.cover <- colonizers*avg.new.cover
#   tmp.out[zeros] <- colonizer.cover
#   
#   ex.mat[t,] <- tmp.out
# }
# rcp85 <- ex.mat[(burn.in+1):time.steps, ]
# 
# 
# ####
# ####  Plot climate change results
# ####
# all.sims <- data.frame(pixel=c(1:length(eq.pixels)),
#                        Baseline=eq.pixels,
#                        RCP45=colMeans(rcp45),
#                        RCP60=colMeans(rcp60),
#                        RCP85=colMeans(rcp85))
# all4plot <- melt(all.sims, id.vars = "pixel")
# 
# 
# 
# 
# ####
# ####  Plot spatial projections for RCP scenarios
# ####
# proj.equil <- data.frame(Lon=subset(growD, Year==1985)$Lon, 
#                          Lat=subset(growD, Year==1985)$Lat,
#                          CURRENT=eq.pixels,
#                          RCP45=colMeans(rcp45),
#                          RCP60=colMeans(rcp60),
#                          RCP85=colMeans(rcp85))
# proj.equil2 <- melt(proj.equil, id.vars = c("Lon", "Lat"))
# 
# # png("../results/climchange_small.png", width = 6, height=4.5, units = "in", res=150)
# ggplot(proj.equil2, aes(x=Lon, y=Lat))+
#   geom_raster(aes(z=value, fill=value))+
#   scale_fill_gradientn(colours=myPalette(200), name="% Cover")+
#   facet_wrap("variable", ncol=1)+
#   coord_equal()+
#   tmp.theme+
#   theme(strip.background=element_rect(fill="white"))
# ggsave("../results/clim_change_mean_spatial.png", height=8, width=4)
# 
# # Make plot of differences to CURRENT
# proj.cast <- dcast(proj.equil2, Lon+Lat~variable)
# projs.only <- proj.cast[,which(colnames(proj.cast) %in% c("RCP45", "RCP60", "RCP85"))]
# diffs.only <- projs.only - proj.cast$CURRENT
# diffs.only$Lon <- proj.cast$Lon
# diffs.only$Lat <- proj.cast$Lat
# diffs.df <- melt(diffs.only, id.vars = c("Lon", "Lat"))
# etaPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
# ggplot(diffs.df, aes(x=Lon, y=Lat))+
#   geom_raster(aes(z=value, fill=value))+
#   scale_fill_gradientn(colours=etaPalette(100), name="% Cover",
#                        limits=c(-7, 7))+
#   facet_wrap("variable", ncol=1)+
#   coord_equal()+
#   tmp.theme+
#   theme(strip.background=element_rect(fill="white"))
# ggsave("../results/clim_change_diffs_spatial.png", height=6, width=4)

# ggplot(all4plot)+
# #   geom_line(stat="density", aes(x=value, y=..density.., fill=variable), size=1, alpha=0.5)+
#   geom_density(aes(x=value, y=..density.., fill=variable, alpha=variable, color=variable), size=1, adjust=2)+
#   scale_fill_manual(values=c("dodgerblue","tan","coral","darkred"), name="Scenario") +
#   scale_color_manual(values=c("dodgerblue","tan","coral","darkred"), name="Scenario") +
#   scale_alpha_manual(values=c(1,0.75,0.5,0.5))+
#   guides(alpha=FALSE)+
#   xlab("Equilibrium Percent Cover")+
#   ylab("Estimated Kernel Density")
# ggsave("../results/clim_change_densities.png", height=4, width=5, dpi = 200)
# gout <- grid.arrange(g1, g2, ncol=1, nrow=2)
# dev.off()







###########################################################
###########################################################
###########################################################
####                                                   ####
####  Run population model with parameter uncertainty  ####
####                                                   ####
###########################################################
###########################################################
###########################################################
# 
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
#       dens.dep <- beta_mu*log(ex.arr[counter,t-1,])
#       tmp.mu <- int_mu + dens.dep + sum(betas*Xtmp)
#       tmp.mu <- exp(tmp.mu + eta)
#       tmp.out <- rpois(ncol(ex.arr[counter,,]), lambda = tmp.mu)
#       
#       #Colonization
#       zeros <- which(ex.arr[counter,t-1,]==0)
#       colonizers <- rbinom(length(zeros), size = 1, antilogit(col.intercept))
#       colonizer.cover <- colonizers*avg.new.cover
#       tmp.out[zeros] <- colonizer.cover
#       
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
