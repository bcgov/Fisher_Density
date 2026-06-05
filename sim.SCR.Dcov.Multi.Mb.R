sim.SCR.Dcov.Multi.Mb <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
           p0=NA,sigma=NA,
           K=NA,X=NA,
           xlim=NA,ylim=NA,res=NA){
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.SCR.Dcov.Mb(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],InSS=InSS[[g]],
                                          p0=p0[g,],sigma=sigma[g],
                                          K=K[g],X=X[[g]],xlim=xlim[g,],ylim=ylim[g,],res=res[g])
    }
    return(data)
  }