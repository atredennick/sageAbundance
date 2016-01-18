n.sim=100000
mu=2.8
sig=0.3

lam=exp(rnorm(n.sim,mu,sig))
y.1=rpois(n.sim,lam)
y.2=rpois(n.sim,exp(mu))

mean(y.1)
mean(y.2)

library(plotrix)
multhist(list(y.1,y.2),breaks=200,col=c(1,2),
         xlab="percent cover", ylab="frequency", border=c(1,2))
text(210, 7000,"approx. posterior", col=2)
text(230, 2000,"true posterior", col=1)
