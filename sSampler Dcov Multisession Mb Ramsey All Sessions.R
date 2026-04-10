#for sessions with cameras
sSamplerCams <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    J1 <- control$J1
    J2 <- control$J2
    K1 <- control$K1
    res <- control$res
    xlim <- control$xlim
    ylim <- control$ylim
    n.cells <- control$n.cells
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    s.nodes <- c(model$expandNodeNames(paste("s[",g,",",i,",",1:2,"]")),
                 model$expandNodeNames(paste("s.cell[",g,",",i,"]")),
                 model$expandNodeNames(paste("dummy.data[",g,",",i,"]")))
    y1.nodes <- model$expandNodeNames(paste("y1[",g,",",i,",",1:J1,",",1:K1,"]"))
    y2.nodes <- model$expandNodeNames(paste("y2[",g,",1:",J2,"]"))
    kern.nodes <- model$expandNodeNames(paste("kern[",g,",",i,",",1:J1,"]"))
    pd1.p.nodes <- model$expandNodeNames(paste("pd1.p[",g,",",i,",",1:J1,"]"))
    pd1.c.nodes <- model$expandNodeNames(paste("pd1.c[",g,",",i,",",1:J1,"]"))
    pd2.nodes <- model$expandNodeNames(paste("pd2[",g,",",i,",",1:J2,"]"))
    pd2.j.nodes <- model$expandNodeNames(paste("pd2.j[",g,"]"))
    # calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    # isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    # calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    # calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
    # if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    # if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    # if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run <- function() {
    z <- model$z[g,i]
    if(z==0){ #propose from prior
      #propose new cell
      model$s.cell[g,i] <<- rcat(1,model$pi.cell[g,1:n.cells])
      #propose x and y in new cell
      s.cell.x <- model$s.cell[g,i]%%n.cells.x 
      s.cell.y <- floor(model$s.cell[g,i]/n.cells.x)+1
      if(s.cell.x==0){
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y-1
      }
      xlim.cell <- c(s.cell.x-1,s.cell.x)*res
      ylim.cell <- c(s.cell.y-1,s.cell.y)*res
      model$s[g,i, 1:2] <<- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
      model$calculate(s.nodes)
      copy(from = model, to = mvSaved, row = 1, nodes = s.nodes, logProb = TRUE)
    }else{ #MH
      s.cand <- c(rnorm(1,model$s[g,i,1],scale), rnorm(1,model$s[g,i,2],scale))
      inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        # s.nodes
        # y1.nodes
        # y2.nodes
        # kern.nodes
        # pd1.p.nodes
        # pd1.c.nodes
        # pd2.nodes
        # pd2.j.nodes
        # browser()
        
        #get initial logprobs
        lp_initial_s <- model$getLogProb(s.nodes)
        lp_initial_y1 <- model$getLogProb(y1.nodes)
        lp_initial_y2 <- model$getLogProb(y2.nodes)
        #pull this out of model object
        nondetect.probs.initial <- 1 - model$pd2.j[g,1:J2] #p(no detect)
        #update proposed s
        model$s[g,i, 1:2] <<- s.cand
        lp_proposed_s <- model$calculate(s.nodes) #proposed logprob for s.nodes
        model$calculate(kern.nodes) #update kern nodes
        model$calculate(pd1.p.nodes) #update pd1.p nodes
        model$calculate(pd1.c.nodes) #update pd1.c nodes
        lp_proposed_y1 <- model$calculate(y1.nodes) #get proposed y1 logProb
        #before updating pd2 nodes, divide out initial 1-pd2[g,i,] from nondetect.probs
        nondetect.probs.proposed <- nondetect.probs.initial/(1-model$pd2[g,i,1:J2])
        model$calculate(pd2.nodes) #update pd2 nodes
        #now multiply nondetect.probs by proposed 1-pd2[g,i,]
        nondetect.probs.proposed <- nondetect.probs.proposed*(1-model$pd2[g,i,1:J2])
        #put these in model object
        model$pd2.j[g,1:J2] <<- 1 - nondetect.probs.proposed
        lp_proposed_y2 <- model$calculate(y2.nodes) #get proposed y2 logProb
        lp_initial <- lp_initial_s + lp_initial_y1 + lp_initial_y2
        lp_proposed <- lp_proposed_s + lp_proposed_y1 + lp_proposed_y2
        log_MH_ratio <- lp_proposed - lp_initial
        accept <- decide(log_MH_ratio)
        if(accept) {
          copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
          copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive){ #we only tune for z=0 proposals
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        # if(reflective) {
        #   lower <- model$getBound(target, 'lower')
        #   upper <- model$getBound(target, 'upper')
        #   if(scale >= 0.5*(upper-lower)) {
        #     scale <<- 0.5*(upper-lower)
        #   }
        # }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)

#for sessions without cams
sSamplerNoCams <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    J1 <- control$J1
    K1 <- control$K1
    res <- control$res
    xlim <- control$xlim
    ylim <- control$ylim
    n.cells <- control$n.cells
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    s.nodes <- c(model$expandNodeNames(paste("s[",g,",",i,",",1:2,"]")),
                 model$expandNodeNames(paste("s.cell[",g,",",i,"]")),
                 model$expandNodeNames(paste("dummy.data[",g,",",i,"]")))
    y1.nodes <- model$expandNodeNames(paste("y1[",g,",",i,",",1:J1,",",1:K1,"]"))
    kern.nodes <- model$expandNodeNames(paste("kern[",g,",",i,",",1:J1,"]"))
    pd1.p.nodes <- model$expandNodeNames(paste("pd1.p[",g,",",i,",",1:J1,"]"))
    pd1.c.nodes <- model$expandNodeNames(paste("pd1.c[",g,",",i,",",1:J1,"]"))
    # calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    # isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    # calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    # calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
    # if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    # if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    # if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run <- function() {
    z <- model$z[g,i]
    if(z==0){ #propose from prior
      #propose new cell
      model$s.cell[g,i] <<- rcat(1,model$pi.cell[g,1:n.cells])
      #propose x and y in new cell
      s.cell.x <- model$s.cell[g,i]%%n.cells.x 
      s.cell.y <- floor(model$s.cell[g,i]/n.cells.x)+1
      if(s.cell.x==0){
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y-1
      }
      xlim.cell <- c(s.cell.x-1,s.cell.x)*res
      ylim.cell <- c(s.cell.y-1,s.cell.y)*res
      model$s[g,i, 1:2] <<- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
      model$calculate(s.nodes)
      copy(from = model, to = mvSaved, row = 1, nodes = s.nodes, logProb = TRUE)
    }else{ #MH
      s.cand <- c(rnorm(1,model$s[g,i,1],scale), rnorm(1,model$s[g,i,2],scale))
      inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        #get initial logprobs
        lp_initial_s <- model$getLogProb(s.nodes)
        lp_initial_y1 <- model$getLogProb(y1.nodes)
        #update proposed s
        model$s[g,i, 1:2] <<- s.cand
        lp_proposed_s <- model$calculate(s.nodes) #proposed logprob for s.nodes
        model$calculate(kern.nodes) #update kern nodes
        model$calculate(pd1.p.nodes) #update pd1.p nodes
        model$calculate(pd1.c.nodes) #update pd1.c nodes
        lp_proposed_y1 <- model$calculate(y1.nodes) #get proposed y1 logProb
        lp_initial <- lp_initial_s + lp_initial_y1
        lp_proposed <- lp_proposed_s + lp_proposed_y1
        log_MH_ratio <- lp_proposed - lp_initial
        accept <- decide(log_MH_ratio)
        if(accept) {
          copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
          copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive){ #we only tune for z=0 proposals
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        # if(reflective) {
        #   lower <- model$getBound(target, 'lower')
        #   upper <- model$getBound(target, 'upper')
        #   if(scale >= 0.5*(upper-lower)) {
        #     scale <<- 0.5*(upper-lower)
        #   }
        # }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)