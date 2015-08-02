model{
  #process model
  for(t in 2:nYears){
    for(i in 1:nCells){
      lambda[i,t] <- exp(intMu + betaMu*lambda[i,t-1] + eta[i])
    }
  }
  
  #initial conditions
  for(k in 1:nCells){
    lambda[k,1] <- exp(intMu + betaMu*lambda0 + eta[k])
  }
  
  #expand alpha effects for full spatial grid
  eta <- K%*%alpha
  
  #likelihood
  for(i in 1:nObs){
    C[i] ~ dpois(lambda[cellMod[i], years[i]])
  }
  
  #alpha priors (iid)
  for(j in 1:nKnots){
    alpha[j] ~ dnorm(0,tau)
  }
  
  #parameter priors
  lambda0 ~ dunif(0,10)
  intMu ~ dnorm(0,0.01)
  betaMu ~ dnorm(0,0.01)
  tau ~ dgamma(1,0.01)
}
