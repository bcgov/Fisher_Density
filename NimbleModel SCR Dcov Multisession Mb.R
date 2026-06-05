NimModel <- nimbleCode({
  sigma.fixed ~ dunif(0,100) #share sigma across sessions
  beta0.p0.SCR.fixed ~ dlogis(0,1) # share first capture probability
  beta1.p0.SCR.fixed ~ dnorm(0,sd=10) #share Mb offset for recapture probability increase from baseline
  for(g in 1:N.session){
    #--------------------------------------------------------------
    # priors
    #--------------------------------------------------------------
    #Density covariates
    D0[g] ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
    # D.beta0 ~ dnorm(0,sd=10)
    D.beta1[g] ~ dnorm(0,sd=10)
    #detection priors
    logit(p0[g,1]) <- beta0.p0.SCR.fixed
    logit(p0[g,2]) <- beta0.p0.SCR.fixed + beta1.p0.SCR.fixed
    sigma[g] <- sigma.fixed
    
    #--------------------------------------------------------------
    #Density model
    D.intercept[g] <- D0[g]*cellArea[g]
    # D.intercept <- exp(D.beta0)*cellArea
    lambda.cell[g,1:n.cells[g]] <- InSS[g,1:n.cells[g]]*exp(D.beta1[g]*D.cov[g,1:n.cells[g]])
    pi.denom[g] <- sum(lambda.cell[g,1:n.cells[g]])
    pi.cell[g,1:n.cells[g]] <- lambda.cell[g,1:n.cells[g]]/pi.denom[g] #expected proportion of total N in cell c
    lambda.N[g] <- D.intercept[g]*pi.denom[g] #Expected N
    N[g] ~ dpois(lambda.N[g]) #realized N in state space
    for(i in 1:M[g]){
      #dunif() here implies uniform distribution within a grid cell
      #also tells nimble s's are in continuous space, not discrete
      s[g,i,1] ~  dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~  dunif(ylim[g,1],ylim[g,2])
      #get cell s_i lives in using look-up table
      s.cell[g,i] <- cells[g,trunc(s[g,i,1]/res[g])+1,trunc(s[g,i,2]/res[g])+1]
      #categorical likelihood for this cell, equivalent to zero's trick
      #also disallowing s's in non-habitat
      dummy.data[g,i] ~ dCell(pi.cell[g,s.cell[g,i]],InSS=InSS[g,s.cell[g,i]])
      #SCR Observation model, skipping z_i=0 calculations
      kern[g,i,1:J[g]] <- GetKern(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], z=z[g,i])
      pd.p[g,i,1:J[g]] <- GetDetectionProb(kern = kern[g,i,1:J[g]], p0=p0[g,1], J=J[g], z=z[g,i])
      pd.c[g,i,1:J[g]] <- GetDetectionProb(kern = kern[g,i,1:J[g]], p0=p0[g,2], J=J[g], z=z[g,i])
      y[g,i,1:J[g],1:K[g]] ~ dBinomialMatrix(pd.p=pd.p[g,i,1:J[g]],pd.c=pd.c[g,i,1:J[g]],
                                        y.state=y.state[g,i,1:J[g],1:K[g]],
                                        K2D=K2D[g,1:J[g],1:K[g]],z=z[g,i]) #vectorized obs mod
    }
  }
})
#custom Metropolis-Hastings update for N/z
