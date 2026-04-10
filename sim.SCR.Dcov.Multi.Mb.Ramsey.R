sim.SCR.Dcov.Multi.Mb.Ramsey <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
           p0.SCR=NA,p0.USCR=NA,sigma=NA,
           K1=NA,K2=NA,X1=NA,X2=NA,
           xlim=NA,ylim=NA,res=NA){
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.SCR.Dcov.Mb.Ramsey(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],InSS=InSS[[g]],
                                          p0.SCR=p0.SCR[g,],p0.USCR=p0.USCR[g],sigma=sigma[g],
                                          K1=K1[g],K2=K2[g],X1=X1[[g]],X2=X2[[g]],
                                          xlim=xlim[g,],ylim=ylim[g,],res=res[g])
    }
    return(data)
  }