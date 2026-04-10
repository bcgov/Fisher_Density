e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.Dcov.Mb.Ramsey <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,lam0=NA,
           p0.SCR=NA,p0.USCR=NA,sigma=NA,
           K1=NA,K2=NA,X1=NA,X2=NA,xlim=NA,ylim=NA,res=NA){
    #get expected N
    cellArea <- res^2
    lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
    lambda.N <- sum(lambda.cell)
    #simulate realized N
    N <- rpois(1,lambda.N)
    
    #recreate some Dcov things so we can pass fewer arguments into this function
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    
    # simulate a population of activity centers
    pi.cell <- lambda.cell/sum(lambda.cell)
    #zero out non-habitat
    pi.cell[InSS==0] <- 0
    s.cell <- sample(1:n.cells,N,prob=pi.cell,replace=TRUE)
    #distribute activity centers uniformly inside cells
    s <- matrix(NA,nrow=N,ncol=2)
    for(i in 1:N){
      tmp <- which(cells==s.cell[i],arr.ind=TRUE) #x and y number
      s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]+res/2])
      s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]+res/2])
    }
    D1 <- e2dist(s,X1)
    J1 <- nrow(X1)
    D2 <- e2dist(s,X2)
    J2 <- nrow(X2)
    
    # Capture individuals - SCR
    y1 <- array(0,dim=c(N,J1,K1))
    y1.state <- array(0,dim=c(N,J1,K1))
    pd1.p <- p0.SCR[1]*exp(-D1*D1/(2*sigma*sigma))
    pd1.c <- p0.SCR[2]*exp(-D1*D1/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J1){
        for(k in 1:K1){
          y1[i,j,k] <- rbinom(1,1,pd1.p[i,j]*(y1.state[i,j,k]==0) + 
                                pd1.c[i,j]*(y1.state[i,j,k]==1))
          if(y1[i,j,k]==1&k!=K1){
            y1.state[i,j,(k+1):K1] <- 1
          }
        }
      }
    }
    
    # Capture individuals - USCR
    y2 <- array(0,dim=c(N,J2,K2))
    pd2 <- p0.USCR*exp(-D2*D2/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J2){
        for(k in 1:K2){
          y2[i,j,k] <- rbinom(1,1,pd2[i,j])
        }
      }
    }
    
    #Get SCR observed data
    caught1 <- which(apply(y1,c(1),sum)>0)
    n1 <- length(caught1)
    # y.SCR.true <- y1
    y1 <- y1[caught1,,]
    if(K1==1){
      y1 <- array(y1,dim=c(dim(y1),1))
    }
    
    #get USCR observed data
    caught2 <- which(apply(y2,c(1),sum)>0)
    n2 <- length(caught2)
    # y.USCR.true <- y2
    y2 <- 1*(apply(y2,c(2,3),sum)>0)
    
    #create trap operation. assumed to be complete in data simulator
    K2D1 <- matrix(1,J1,K1)
    K1D2 <- rep(K2,J2)
    
    out <- list(y1=y1,y2=y2,X1=X1,X2=X2,K1=K1,K2=K2,
                K2D1=K2D1,K1D2=K1D2,s=s,n1=n1,n2=n2,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N)
    return(out)
  }
