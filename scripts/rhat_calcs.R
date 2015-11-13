##  Script to read in MCMC results and calculate
##  scale reduction factors

rm(list=ls())

####
####  Load libraries -----------------------------------------------------------
####
library(ggplot2)
library(dplyr)


####
####  Define Rhat function; stolen from `ggmcmc` ggs_Rhats() fxn ---------------
####
get_rhats <- function (D, family = NA, scaling = 1.5) 
{
  if (attributes(D)$nChains < 2) {
    stop("At least two chains are required")
  }
  if (!is.na(family)) {
    D <- get_family(D, family = family)
  }
  psi.dot <- D %>% group_by(Parameter, Chain) %>% dplyr::summarize(psi.dot = mean(value))
  psi.j <- D %>% group_by(Parameter) %>% dplyr::summarize(psi.j = mean(value))
  b.df <- dplyr::inner_join(psi.dot, psi.j, by = "Parameter")
  B <- b.df %>% group_by(Parameter) %>% dplyr::summarize(B = var(psi.j - 
                                                                   psi.dot) * attributes(D)$nIterations)
  B <- unique(B)
  s2j <- D %>% group_by(Parameter, Chain) %>% dplyr::summarize(s2j = var(value))
  W <- s2j %>% group_by(Parameter) %>% dplyr::summarize(W = mean(s2j))
  BW <- dplyr::inner_join(B, W, by = "Parameter") %>% dplyr::mutate(wa = (((attributes(D)$nIterations - 
                                                                              1)/attributes(D)$nIterations) * W) + ((1/attributes(D)$nIterations) * 
                                                                                                                      B), Rhat = sqrt(wa/W))
  BW$Rhat[is.nan(BW$Rhat)] <- NA
  f <- ggplot(BW, aes(x = Rhat, y = Parameter)) + geom_point() + 
    xlab(expression(hat("R"))) + ggtitle("Potential Scale Reduction Factor")
  if (!is.na(scaling)) {
    scaling <- ifelse(scaling > max(BW$Rhat, na.rm = TRUE), 
                      scaling, max(BW$Rhat, na.rm = TRUE))
    f <- f + xlim(min(BW$Rhat, na.rm = TRUE), scaling)
  }
  return(BW)
}


####
####  Read in MCMC 
####
setwd("/Users/atredenn/Repos/sageAbundance/results/")
outmcmcfile <- "poissonSage_randYear_mcmc.RDS"
out.mcmc <- readRDS(outmcmcfile)
colnames(out.mcmc) <- c("Iteration", "Parameter", "value", "mcmc_iter", "Chain")
attributes(out.mcmc)$nChains <- length(unique(out.mcmc$Chain))
attributes(out.mcmc)$nIterations <- length(unique(out.mcmc$Iteration))
rhats <- get_rhats(out.mcmc)
rhats$facet_id <- rep(c(1,2), each=nrow(rhats)/2)

ggplot(rhats, aes(x=Rhat, y=Parameter))+
  geom_point()+
  geom_vline(aes(xintercept=1.1), col="red")+
  facet_wrap("facet_id", ncol=2, scales="free_y")+
  xlim(min(rhats$Rhat, na.rm = TRUE), 1.5)+
  xlab(expression(hat("R"))) + 
  ggtitle("Potential Scale Reduction Factors")
  
ggsave("rhat_plot.png", height = 10, dpi = 100)


