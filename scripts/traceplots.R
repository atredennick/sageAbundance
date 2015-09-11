
setwd("/Users/atredenn/Desktop/")
maxiter <- 5
longchain1 <- readRDS("POIS_iterchunk1_chain2.RDS")
for(i in 2:maxiter){
  longchain1 <- rbind(longchain1, readRDS(paste0("POIS_iterchunk",i,"_chain2.RDS")))
}

longchain2 <- readRDS("POIS_iterchunk1_chain3.RDS")
for(i in 2:maxiter){
  longchain2 <- rbind(longchain2, readRDS(paste0("POIS_iterchunk",i,"_chain3.RDS")))
}

longchain3 <- readRDS("POIS_iterchunk1_chain1.RDS")
for(i in 2:maxiter){
  longchain3 <- rbind(longchain3, readRDS(paste0("POIS_iterchunk",i,"_chain1.RDS")))
}


all.params <- unique(longchain1$Parameter)
iters <- c(100:1000)
iters <- iters[-which(iters %in% c(201:205,401:405,601:605,801:805))]

pdf("test.pdf", width = 8.5, height = 11)
par(mfrow=c(6,3))
for(plot.param in all.params){
  plot(unlist(longchain1[which(longchain1$Parameter==plot.param),"value"])[iters], 
       type="l", ylab="value", main=plot.param, xlab="iteration")
  lines(unlist(longchain2[which(longchain2$Parameter==plot.param),"value"])[iters], col="darkorange")
  lines(unlist(longchain3[which(longchain3$Parameter==plot.param),"value"])[iters], col="purple")
}
dev.off()


