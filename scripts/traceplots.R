##  Script to read in MCMC results from the DoRC server
##  Creates traceplots and save MCMC output for inference (1000 iterations)

#Poisson fit
filestart <- "POIS_iterchunk" 
traceplotfile <- "poissonSage_randYear_traceplots.pdf"
outmcmcfile <- "poissonSage_randYear_mcmc.RDS"

#Negative binomial fit
# filestart <- "iterchunk" 
# traceplotfile <- "negbinomSage_traceplots.pdf"
# outmcmcfile <- "negbinomSage_mcmc.RDS"

####
####  Load libraries
####
library(ggplot2)
library(scales)


####
####  Read in MCMC chunks from server
####
setwd("/Volumes/A02046115/sage2/")
maxiter <- 10
longchain1 <- readRDS(paste0(filestart, "1_chain2.RDS"))
for(i in 2:maxiter){
  longchain1 <- rbind(longchain1, readRDS(paste0(filestart,i,"_chain2.RDS")))
}

longchain2 <- readRDS(paste0(filestart, "1_chain3.RDS"))
for(i in 2:maxiter){
  longchain2 <- rbind(longchain2, readRDS(paste0(filestart,i,"_chain3.RDS")))
}

longchain3 <- readRDS(paste0(filestart, "1_chain1.RDS"))
for(i in 2:maxiter){
  longchain3 <- rbind(longchain3, readRDS(paste0(filestart,i,"_chain1.RDS")))
}



####
####  Create subset vectors to get 1000 posterior samples
####
iters <- c(976:2000) # equals 1000 samples after next line
## throw away 5 iterations from the beginning of each chunk (warmup in STAN)
iters <- iters[-which(iters %in% c(1001:1005,1201:1205,1401:1405,
                                   1601:1605,1801:1805))]


####
####  Save plot of traceplots for each parameter
####
setwd("/Users/atredenn/Repos/sageAbundance/results/")

all.params <- unique(longchain1$Parameter)

pdf(traceplotfile, width = 8.5, height = 11)
par(mfrow=c(6,3))
for(plot.param in all.params){
  plot(unlist(longchain1[which(longchain1$Parameter==plot.param),"value"])[iters], 
       type="l", ylab="value", main=plot.param, xlab="iteration", col=alpha("black",0.85))
  lines(unlist(longchain2[which(longchain2$Parameter==plot.param),"value"])[iters], col=alpha("darkorange",0.75))
  lines(unlist(longchain3[which(longchain3$Parameter==plot.param),"value"])[iters], col=alpha("purple",0.65))
}
dev.off()



####
####  Make one long MCMC dataframe for inference and modeling
####
lc1 <- longchain1[with(longchain1, order(Parameter)), ]
lc1$mcmc_iter <- rep(c(1:2000), length(unique(lc1$Parameter)))
lc1$chain <- 1
lc1 <- lc1[which(lc1$mcmc_iter %in% iters),]
lc1$mcmc_iter <- rep(c(1:1000), length(unique(lc1$Parameter)))

lc2 <- longchain2[with(longchain2, order(Parameter)), ]
lc2$mcmc_iter <- rep(c(1:2000), length(unique(lc2$Parameter)))
lc2$chain <- 2
lc2 <- lc2[which(lc2$mcmc_iter %in% iters),]
lc2$mcmc_iter <- rep(c(1:1000), length(unique(lc2$Parameter)))

lc3 <- longchain3[with(longchain2, order(Parameter)), ]
lc3$mcmc_iter <- rep(c(1:2000), length(unique(lc3$Parameter)))
lc3$chain <- 3
lc3 <- lc3[which(lc3$mcmc_iter %in% iters),]
lc3$mcmc_iter <- rep(c(1:1000), length(unique(lc3$Parameter)))

lc.save <- rbind(lc1, lc2, lc3)

## Just for testingm, commented out now
# ggplot(subset(lc.save, Parameter=="beta_mu"), aes(x=mcmc_iter, y=value, color=as.factor(chain)))+
#   geom_line()
# max(lc.save$mcmc_iter)
# min(lc.save$mcmc_iter)



####
####  Save the MCMC for inference
####
setwd("/Users/atredenn/Repos/sageAbundance/results/")
saveRDS(lc.save, outmcmcfile)


