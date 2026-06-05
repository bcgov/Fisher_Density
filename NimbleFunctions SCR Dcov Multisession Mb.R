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

GetDetectionProb <- nimbleFunction(
  run = function(kern=double(1),p0=double(0),J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      ans <- p0*kern
      return(ans)
    }
  }
)

dBinomialMatrix <- nimbleFunction(
  run = function(x = double(2), pd.p = double(1), pd.c=double(1), y.state = double(2), K2D = double(2), z = double(0),
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
            if(y.state[j,k]==0){
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
  run = function(n = integer(0),pd.p = double(1), pd.c=double(1), y.state = double(2), K2D=double(2), z = double(0)) {
    returnType(double(2))
    J <- nimDim(pd.p)[1]
    K <- nimDim(pd.p)[2]
    out <- matrix(nrow=J,ncol=K,value=0)
    return(out)
  }
)

#Required custom update for number of calls
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("N","y","kern","pd.p","pd.c")) #nodes to copy back to mvSaved
    inds.detected <- control$inds.detected
    z.ups <- control$z.ups
    J <- control$J
    K <- control$K
    M <- control$M
    g <- control$g
    #nodes used for update, calcNodes + z nodes
    y.nodes <- model$expandNodeNames(paste("y[",g,",",1:M,",",1:J,",",1:K,"]"))
    kern.nodes <- model$expandNodeNames(paste("kern[",g,",",1:M,",",1:J,"]"))
    pd.p.nodes <- model$expandNodeNames(paste("pd.p[",g,",",1:M,",",1:J,"]"))
    pd.c.nodes <- model$expandNodeNames(paste("pd.c[",g,",",1:M,",",1:J,"]"))
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
          lp.initial.y <- model$getLogProb(y.nodes[pick])

          #propose new N/z
          model$N[g] <<-  model$N[g] - 1
          model$z[g,pick] <<- 0

          #turn pd off
          #don't use calculate for pd2.j, compute and insert manually
          model$calculate(kern.nodes[pick])
          model$calculate(pd.p.nodes[pick])
          model$calculate(pd.c.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            mvSaved["kern",1][g,pick,1:J] <<- model[["kern"]][g,pick,1:J]
            mvSaved["pd.p",1][g,pick,1:J] <<- model[["pd.p"]][g,pick,1:J]
            mvSaved["pd.c",1][g,pick,1:J] <<- model[["pd.c"]][g,pick,1:J]
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            model[["kern"]][g,pick,1:J] <<- mvSaved["kern",1][g,pick,1:J]
            model[["pd.p"]][g,pick,1:J] <<- mvSaved["pd.p",1][g,pick,1:J]
            model[["pd.c"]][g,pick,1:J] <<- mvSaved["pd.c",1][g,pick,1:J]
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y.nodes[pick])
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
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[g] <<-  model$N[g] + 1
          model$z[g,pick] <<- 1
          
          #turn pd on
          model$calculate(kern.nodes[pick])
          model$calculate(pd.p.nodes[pick])
          model$calculate(pd.c.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            mvSaved["kern",1][g,pick,1:J] <<- model[["kern"]][g,pick,1:J]
            mvSaved["pd.p",1][g,pick,1:J] <<- model[["pd.p"]][g,pick,1:J]
            mvSaved["pd.c",1][g,pick,1:J] <<- model[["pd.c"]][g,pick,1:J]
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            model[["kern"]][g,pick,1:J] <<- mvSaved["kern",1][g,pick,1:J]
            model[["pd.p"]][g,pick,1:J] <<- mvSaved["pd.p",1][g,pick,1:J]
            model[["pd.c"]][g,pick,1:J] <<- mvSaved["pd.c",1][g,pick,1:J]
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y.nodes[pick])
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