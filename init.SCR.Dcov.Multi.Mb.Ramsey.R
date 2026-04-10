init.SCR.Dcov.Multi.Mb.Ramsey <-
  function(data=data,M=M){
    N.session <- length(data)
    if(length(M)!=N.session)stop("M and data must be of length 'N.session'")
    M.max <- max(M)
    J1 <- unlist(lapply(data,function(x){nrow(x$X1)}))
    J2.tmp <- lapply(data,function(x){nrow(x$X2)})
    has.cams <- which(unlist(lapply(J2.tmp,is.null)))
    if(length(has.cams)>0){
      J2.tmp[has.cams] <- 0
    }
    J2 <- unlist(J2.tmp)
    J1.max <- max(J1)
    J2.max <- max(J2)
    K1 <- unlist(lapply(data,function(x){x$K1}))
    K2 <- unlist(lapply(data,function(x){x$K2}))
    K1.max <- max(K1)
    K2.max <- max(K2)
    n.cells <- unlist(lapply(data,function(x){x$n.cells}))
    n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
    n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
    n.cells.max <- max(n.cells)
    n.cells.x.max <- max(n.cells.x)
    n.cells.y.max <- max(n.cells.y)
    #Structure data for nimble
    y1 <- array(0,dim=c(N.session,M.max,J1.max,K1.max)) #maximal augmentation across sessions
    y2 <- array(0,dim=c(N.session,J2.max))
    X1 <- array(0,dim=c(N.session,J1.max,2))
    X2 <- array(0,dim=c(N.session,J2.max,2))
    xlim <- ylim <- matrix(0,N.session,2)
    n1 <- rep(NA,N.session)
    K2D1 <- array(0,dim=c(N.session,J1.max,K1.max))
    K1D2 <- array(0,dim=c(N.session,J2.max))
    
    res <- unlist(lapply(data,function(x){x$res}))
    cellArea <- res^2
    x.vals <- matrix(NA,N.session,n.cells.x.max)
    y.vals <- matrix(NA,N.session,n.cells.y.max)
    dSS <- array(NA,dim=c(N.session,n.cells.max,2))
    InSS <- array(0,dim=c(N.session,n.cells.max))
    D.cov <- array(NA,dim=c(N.session,n.cells.max))
    cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
    for(g in 1:N.session){
      y1[g,1:data[[g]]$n1,1:J1[g],1:K1[g]] <- data[[g]]$y1
      X1[g,1:J1[g],1:2] <- data[[g]]$X1
      if(!is.na(data[[g]]$K2)){
        y2[g,1:J2[g]] <- rowSums(data[[g]]$y2)
        X2[g,1:J2[g],1:2] <- data[[g]]$X2
        K1D2[g,1:J2[g]] <- data[[g]]$K1D2
      }else{
        print(paste("No camera data in session",g))
      }
      xlim[g,] <- data[[g]]$xlim
      ylim[g,] <- data[[g]]$ylim
      n1[g] <- data[[g]]$n1
      K2D1[g,1:J1[g],1:K1[g]] <- data[[g]]$K2D1
      x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
      y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
      dSS[g,1:n.cells[g],] <- data[[g]]$dSS
      InSS[g,1:n.cells[g]] <- data[[g]]$InSS
      D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
      cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
    }
    
    #get state matrix
    y1.state <- array(0,dim=c(g,M.max,J1.max,K1.max))
    for(g in 1:N.session){
      for(i in 1:M[g]){
        for(j in 1:J1[g]){
          if(sum(y1[g,i,j,])>0){
            firstcap <- which(y1[g,i,j,]>0)[1]
            if(firstcap<K1[g]){
              y1.state[g,i,j,(firstcap+1):K1[g]] <- 1
            }
          }
        }
      }
    }
    
    s.init <- array(NA,dim=c(N.session,M.max,2))
    y1.3D <- apply(y1,c(1,2,3),sum)
    y1.2D <- apply(y1,c(1,2),sum)
    for(g in 1:N.session){
      s.init[g,1:M[g],] <- cbind(runif(M[g],xlim[g,1],xlim[g,2]), runif(M[g],ylim[g,1],ylim[g,2])) #assign random locations
      idx <- which(y1.2D[g,]>0) #switch for those actually caught
      for(i in idx){
        trps <- matrix(X1[g,y1.3D[g,i,]>0,1:2],ncol=2,byrow=FALSE)
        if(nrow(trps)>1){
          s.init[g,i,] <- c(mean(trps[,1]),mean(trps[,2]))
        }else{
          s.init[g,i,] <- trps
        }
      }
    }
    #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
    e2dist  <-  function (x, y){
      i <- sort(rep(1:nrow(y), nrow(x)))
      dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
      matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
    }
    getCell  <-  function(s,res,cells){
      cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
    }
    for(g in 1:N.session){
      alldists <- e2dist(s.init[g,,],data[[g]]$dSS)
      alldists[,data[[g]]$InSS==0] <- Inf
      for(i in 1:M[g]){
        this.cell <- data[[g]]$cells[trunc(s.init[g,i,1]/data[[g]]$res)+1,trunc(s.init[g,i,2]/data[[g]]$res)+1]
        if(data[[g]]$InSS[this.cell]==0){
          cands <- alldists[i,]
          new.cell <- which(alldists[i,]==min(alldists[i,]))
          s.init[g,i,] <- data[[g]]$dSS[new.cell,]
        }
      }
    }
    
    z.init <- 1*(apply(y1,c(1,2),sum)>0)
    #might need to turn on some z's near USCR detectors. T
    #Turn on closest individual to each USCR detector with detections if not already on
    for(g in 1:N.session){
      if(J2[g]>0){ #if session has cameras
        for(j in 1:J2[g]){
          if(y2[g,j]>0){
            dists <- sqrt((X2[g,j,1]-s.init[g,,1])^2 + (X2[g,j,2]-s.init[g,,2])^2)
            idx <- which(dists==min(dists,na.rm=TRUE))
            z.init[g,idx] <- 1
          }
        }
      }
    }
    dummy.data <- matrix(0,N.session,M.max) #dummy data not used, doesn't really matter what the values are
    
    z.data <- matrix(NA,N.session,M.max)
    for(g in 1:N.session){
      z.data[g,1:n1[g]] <- 1
    }
    
    return(list(y1=y1,y1.state=y1.state,y2=y2,X1=X1,X2=X2,xlim=xlim,ylim=ylim,K1=K1,K2=K2,J1=J1,J2=J2,
                n1=n1,K2D1=K2D1,K1D2=K1D2,res=res,cellArea=cellArea,x.vals=x.vals,
                y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
                n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data,
                s.init=s.init,z.init=z.init,z.data=z.data))
  }
