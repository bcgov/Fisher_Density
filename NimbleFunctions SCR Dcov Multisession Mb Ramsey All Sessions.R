dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0), InSS = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(InSS==1){
      logProb <- log(pi.cell)
    }else{
      logProb <- -Inf
    }
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0), InSS = double(0)) {
    returnType(double(0))
    return(0)
  }
)

GetKern <- nimbleFunction(
  run = function(s = double(1), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

GetSCRDetectionProb <- nimbleFunction(
  run = function(kern=double(1),p0=double(0),J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      ans <- p0*kern
      return(ans)
    }
  }
)

GetUSCRDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

GetTrapProbs <- nimbleFunction(
  run = function(pd2 = double(2),z = double(1)){ 
    returnType(double(1))
    J <- nimDim(pd2)[2]
    pd2.j <- rep(0,J)
    z1.idx <- which(z==1) #only consider z=1 inds
    for(j in 1:J){
      pd2.j[j] <- 1 - (prod(1-pd2[z1.idx,j]))
    }
    return(pd2.j)
  }
)

dBinomialMatrix <- nimbleFunction(
  run = function(x = double(2), pd.p = double(1), pd.c=double(1), y1.state = double(2), K2D = double(2), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      J <- nimDim(K2D)[1]
      K <- nimDim(K2D)[2]
      logProb <- matrix(0,J,K)
      for(j in 1:J){
        for(k in 1:K){
          if(K2D[j,k]>0){
            if(y1.state[j,k]==0){
              logProb[j,k] <- dbinom(x[j,k], size = 1, p = pd.p[j], log = TRUE)
            }else{
              logProb[j,k] <- dbinom(x[j,k], size = 1, p = pd.c[j], log = TRUE)
            }
          }
        }
      }
      return(sum(logProb))
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialMatrix <- nimbleFunction(
  run = function(n = integer(0),pd.p = double(1), pd.c=double(1), y1.state = double(2), K2D=double(2), z = double(0)) {
    returnType(double(2))
    J <- nimDim(pd.p)[1]
    K <- nimDim(pd.p)[2]
    out <- matrix(nrow=J,ncol=K,value=0)
    return(out)
  }
)

#Required custom update for number of calls
#this one for camera sessions
zSamplerCams <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("N","y1","y2","kern","pd1.p","pd1.c","pd2","pd2.j")) #nodes to copy back to mvSaved
    inds.detected <- control$inds.detected
    z.ups <- control$z.ups
    J1 <- control$J1
    J2 <- control$J2
    K1 <- control$K1
    M <- control$M
    g <- control$g
    #nodes used for update, calcNodes + z nodes
    y1.nodes <- model$expandNodeNames(paste("y1[",g,",",1:M,",",1:J1,",",1:K1,"]"))
    y2.nodes <- model$expandNodeNames(paste("y2[",g,",1:",J2,"]"))
    kern.nodes <- model$expandNodeNames(paste("kern[",g,",",1:M,",",1:J1,"]"))
    pd1.p.nodes <- model$expandNodeNames(paste("pd1.p[",g,",",1:M,",",1:J1,"]"))
    pd1.c.nodes <- model$expandNodeNames(paste("pd1.c[",g,",",1:M,",",1:J1,"]"))
    pd2.nodes <- model$expandNodeNames(paste("pd2[",g,",",1:M,",",1:J2,"]"))
    pd2.j.nodes <- model$expandNodeNames(paste("pd2.j[",g,"]"))
    N.node <- model$expandNodeNames(paste("N[",g,"]"))
    z.nodes <- model$expandNodeNames(paste("z[",g,",1:",M,"]"))
  },
  run = function() {
    #track these "manually" so computations faster than nimble will do them
    nondetect.probs.initial <- 1 - model$pd2.j[g,1:J2]  #p(no detect)
    # browser()
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z[g,1:M]==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick])
          lp.initial.y2 <- model$getLogProb(y2.nodes)

          #propose new N/z
          model$N[g] <<-  model$N[g] - 1
          model$z[g,pick] <<- 0

          #turn pd off
          #don't use calculate for pd2.j, compute and insert manually
          model$calculate(kern.nodes[pick])
          model$calculate(pd1.p.nodes[pick])
          model$calculate(pd1.c.nodes[pick])
          nondetect.probs.proposed <- nondetect.probs.initial/(1-model$pd2[g,pick,1:J2]) #divide these out before calculate, which sets to 0
          model$calculate(pd2.nodes[pick])
          # model$calculate(pd2.j.nodes)
          model$pd2.j[g,1:J2]  <<- 1 - nondetect.probs.proposed
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick]) #will always be 0
          lp.proposed.y2 <- model$calculate(y2.nodes)

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y1 + lp.proposed.y2) - (lp.initial.N + lp.initial.y1 + lp.initial.y2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            mvSaved["kern",1][g,pick,1:J1] <<- model[["kern"]][g,pick,1:J1]
            mvSaved["pd1.p",1][g,pick,1:J1] <<- model[["pd1.p"]][g,pick,1:J1]
            mvSaved["pd1.c",1][g,pick,1:J1] <<- model[["pd1.c"]][g,pick,1:J1]
            mvSaved["pd2",1][g,pick,1:J2] <<- model[["pd2"]][g,pick,1:J2]
            mvSaved["pd2.j",1][g,1:J2] <<- model[["pd2.j"]][g,1:J2]
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
            nondetect.probs.initial  <- nondetect.probs.proposed
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            model[["kern"]][g,pick,1:J1] <<- mvSaved["kern",1][g,pick,1:J1]
            model[["pd1.p"]][g,pick,1:J1] <<- mvSaved["pd1.p",1][g,pick,1:J1]
            model[["pd1.c"]][g,pick,1:J1] <<- mvSaved["pd1.c",1][g,pick,1:J1]
            model[["pd2"]][g,pick,1:J2] <<- mvSaved["pd2",1][g,pick,1:J2]
            model[["pd2.j"]][g,1:J2] <<- mvSaved["pd2.j",1][g,1:J2]
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y1.nodes[pick])
            model$calculate(y2.nodes)
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[g] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z[g,1:M]==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick]) #will always be 0
          lp.initial.y2 <- model$getLogProb(y2.nodes)
          
          #propose new N/z
          model$N[g] <<-  model$N[g] + 1
          model$z[g,pick] <<- 1
          
          #turn pd on
          model$calculate(kern.nodes[pick])
          model$calculate(pd1.p.nodes[pick])
          model$calculate(pd1.c.nodes[pick])
          model$calculate(pd2.nodes[pick])
          #don't use calculate, compute and insert manually
          # model$calculate(pd2.j.nodes)
          nondetect.probs.proposed <- nondetect.probs.initial*(1-model$pd2[g,pick,1:J2])
          model$pd2.j[g,1:J2]  <<- 1 - nondetect.probs.proposed
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick])
          lp.proposed.y2 <- model$calculate(y2.nodes)
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y1 + lp.proposed.y2) - (lp.initial.N + lp.initial.y1 + lp.initial.y2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            mvSaved["kern",1][g,pick,1:J1] <<- model[["kern"]][g,pick,1:J1]
            mvSaved["pd1.p",1][g,pick,1:J1] <<- model[["pd1.p"]][g,pick,1:J1]
            mvSaved["pd1.c",1][g,pick,1:J1] <<- model[["pd1.c"]][g,pick,1:J1]
            mvSaved["pd2",1][g,pick,1:J2] <<- model[["pd2"]][g,pick,1:J2]
            mvSaved["pd2.j",1][g,1:J2] <<- model[["pd2.j"]][g,1:J2]
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
            nondetect.probs.initial  <- nondetect.probs.proposed
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            model[["kern"]][g,pick,1:J1] <<- mvSaved["kern",1][g,pick,1:J1]
            model[["pd1.p"]][g,pick,1:J1] <<- mvSaved["pd1.p",1][g,pick,1:J1]
            model[["pd1.c"]][g,pick,1:J1] <<- mvSaved["pd1.c",1][g,pick,1:J1]
            model[["pd2"]][g,pick,1:J2] <<- mvSaved["pd2",1][g,pick,1:J2]
            model[["pd2.j"]][g,1:J2] <<- mvSaved["pd2.j",1][g,1:J2]
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y1.nodes[pick])
            model$calculate(y2.nodes)
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # copy(from = model, to = mvSaved, row = 1, nodes = z.nodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#this one for sessions with no cameras
zSamplerNoCams <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("N","y1","kern","pd1.p","pd1.c")) #nodes to copy back to mvSaved
    inds.detected <- control$inds.detected
    z.ups <- control$z.ups
    J1 <- control$J1
    K1 <- control$K1
    M <- control$M
    g <- control$g
    #nodes used for update, calcNodes + z nodes
    y1.nodes <- model$expandNodeNames(paste("y1[",g,",",1:M,",",1:J1,",",1:K1,"]"))
    kern.nodes <- model$expandNodeNames(paste("kern[",g,",",1:M,",",1:J1,"]"))
    pd1.p.nodes <- model$expandNodeNames(paste("pd1.p[",g,",",1:M,",",1:J1,"]"))
    pd1.c.nodes <- model$expandNodeNames(paste("pd1.c[",g,",",1:M,",",1:J1,"]"))
    N.node <- model$expandNodeNames(paste("N[",g,"]"))
    z.nodes <- model$expandNodeNames(paste("z[",g,",1:",M,"]"))
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z[g,1:M]==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick])
          
          #propose new N/z
          model$N[g] <<-  model$N[g] - 1
          model$z[g,pick] <<- 0
          
          #turn pd off
          #don't use calculate for pd2.j, compute and insert manually
          model$calculate(kern.nodes[pick])
          model$calculate(pd1.p.nodes[pick])
          model$calculate(pd1.c.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y1) - (lp.initial.N + lp.initial.y1)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            mvSaved["kern",1][g,pick,1:J1] <<- model[["kern"]][g,pick,1:J1]
            mvSaved["pd1.p",1][g,pick,1:J1] <<- model[["pd1.p"]][g,pick,1:J1]
            mvSaved["pd1.c",1][g,pick,1:J1] <<- model[["pd1.c"]][g,pick,1:J1]
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            model[["kern"]][g,pick,1:J1] <<- mvSaved["kern",1][g,pick,1:J1]
            model[["pd1.p"]][g,pick,1:J1] <<- mvSaved["pd1.p",1][g,pick,1:J1]
            model[["pd1.c"]][g,pick,1:J1] <<- mvSaved["pd1.c",1][g,pick,1:J1]
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y1.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[g] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z[g,1:M]==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick]) #will always be 0

          #propose new N/z
          model$N[g] <<-  model$N[g] + 1
          model$z[g,pick] <<- 1
          
          #turn pd on
          model$calculate(kern.nodes[pick])
          model$calculate(pd1.p.nodes[pick])
          model$calculate(pd1.c.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick])

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y1) - (lp.initial.N + lp.initial.y1)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            mvSaved["kern",1][g,pick,1:J1] <<- model[["kern"]][g,pick,1:J1]
            mvSaved["pd1.p",1][g,pick,1:J1] <<- model[["pd1.p"]][g,pick,1:J1]
            mvSaved["pd1.c",1][g,pick,1:J1] <<- model[["pd1.c"]][g,pick,1:J1]
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            model[["kern"]][g,pick,1:J1] <<- mvSaved["kern",1][g,pick,1:J1]
            model[["pd1.p"]][g,pick,1:J1] <<- mvSaved["pd1.p",1][g,pick,1:J1]
            model[["pd1.c"]][g,pick,1:J1] <<- mvSaved["pd1.c",1][g,pick,1:J1]
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y1.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # copy(from = model, to = mvSaved, row = 1, nodes = z.nodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)