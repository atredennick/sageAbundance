model{
  for(i in 1:nObs){
    C[i] ~ dpois(lambda[i])
    lambda[i] <- exp(intMu + betaMu*N[i] + eta[cell[i]])
  }
  eta <- K%*%alpha
  for(j in 1:nKnots){
    alpha[j] ~ dnorm(0,tau)
  }
  intMu ~ dnorm(0,1e-6)
  betaMu ~ dnorm(0,1e-6)
  tau ~ dgamma(1,0.01)
}
