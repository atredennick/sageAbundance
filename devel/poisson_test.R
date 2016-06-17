n <- 100
mu <- 3
sd.mu <- 0.5
K <- 10000

diff.yhat <- numeric(K)
yhat1avg <- numeric(K)
yhat2avg <- numeric(K)
for(k in 1:K){
  lambda <- exp(rnorm(1, mu, sd.mu))
  y.hat1 <- rpois(n, lambda)
  y.hat2 <- rpois(n, exp(mu))
  diff.yhat[k] <- mean(y.hat1)-mean(y.hat2)
  yhat1avg[k] <- mean(y.hat1)
  yhat2avg[k] <- mean(y.hat2)
}

par(mfrow=c(1,2))
hist(diff.yhat)
plot(yhat1avg, yhat2avg, pch=".")
abline(0,1)


