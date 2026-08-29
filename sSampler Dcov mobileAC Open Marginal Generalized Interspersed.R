#sSampler 1 updates s when z.super=1 and z=1
sSampler1 <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    K.sight <- control$K.sight
    n.marked.all <- control$n.marked.all
    mark.states <- control$mark.states
    mark.states.all <- control$mark.states.all
    M <- dim(mark.states.all)[1]
    s.nodes <- control$s.nodes
    calcNodes <- control$calcNodes
    
    #Explicit observation-model node sets for this primary occasion.
    #Secondary occasions are separated by the focal individual's fixed mark state.
    pd.nodes <- y.mark.nodes <- character(0)
    lam.nodes <- y.mID.nodes <- character(0)
    lam.mnoID.nodes <- y.mnoID.nodes <- character(0)
    lam.um.nodes <- y.um.nodes <- character(0)
    lam.unk.nodes <- y.unk.nodes <- character(0)
    k.marked <- k.unmarked <- integer(0)
    
    if(J.mark[g]>0){
      pd.nodes <- model$expandNodeNames(paste0("pd[",i,",",g,",1:",J.mark[g],"]"))
      y.mark.nodes <- model$expandNodeNames(paste0("y.mark[",i,",",g,",1:",J.mark[g],"]"))
    }
    if(J.sight[g]>0){
      lam.nodes <- model$expandNodeNames(paste0("lam[",i,",",g,",1:",J.sight[g],"]"))
      if(K.sight[g]>0){
        for(k.setup in 1:K.sight[g]){
          if(mark.states[g,k.setup]==1){
            k.marked <- c(k.marked,k.setup)
            if(i<=n.marked.all){
              y.mID.nodes <- c(y.mID.nodes,
                               model$expandNodeNames(paste0("y.mID[",i,",",g,",1:",J.sight[g],",",k.setup,"]")))
            }
            lam.mnoID.nodes <- c(lam.mnoID.nodes,
                                 model$expandNodeNames(paste0("lam.mnoID[",g,",1:",J.sight[g],",",k.setup,"]")))
            y.mnoID.nodes <- c(y.mnoID.nodes,
                               model$expandNodeNames(paste0("y.mnoID[",g,",1:",J.sight[g],",",k.setup,"]")))
          }else{
            k.unmarked <- c(k.unmarked,k.setup)
            lam.um.nodes <- c(lam.um.nodes,
                              model$expandNodeNames(paste0("lam.um[",g,",1:",J.sight[g],",",k.setup,"]")))
            y.um.nodes <- c(y.um.nodes,
                            model$expandNodeNames(paste0("y.um[",g,",1:",J.sight[g],",",k.setup,"]")))
          }
          #unknown-mark-status sightings always depend on the focal individual through
          #either bigLam.marked or bigLam.unmarked
          lam.unk.nodes <- c(lam.unk.nodes,
                             model$expandNodeNames(paste0("lam.unk[",g,",1:",J.sight[g],",",k.setup,"]")))
          y.unk.nodes <- c(y.unk.nodes,
                           model$expandNodeNames(paste0("y.unk[",g,",1:",J.sight[g],",",k.setup,"]")))
        }
      }
    }
    n.pd <- length(pd.nodes)
    n.y.mark <- length(y.mark.nodes)
    n.lam <- length(lam.nodes)
    n.y.mID <- length(y.mID.nodes)
    n.lam.mnoID <- length(lam.mnoID.nodes)
    n.y.mnoID <- length(y.mnoID.nodes)
    n.lam.um <- length(lam.um.nodes)
    n.y.um <- length(y.um.nodes)
    n.lam.unk <- length(lam.unk.nodes)
    n.y.unk <- length(y.unk.nodes)
    n.k.marked <- length(k.marked)
    n.k.unmarked <- length(k.unmarked)
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
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
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function(){
    z.super <- model$z.super[i]
    z <- model$z[i,g]
    if(z.super==1&z==1){
      s.cand <- c(rnorm(1,model$s[i,g,1],scale),rnorm(1,model$s[i,g,2],scale))
      inbox <- s.cand[1]<xlim[2] & s.cand[1]>xlim[1] & s.cand[2]<ylim[2] & s.cand[2]>ylim[1]
      if(inbox){
        #initial log probability for movement/telemetry terms and affected observation likelihoods
        lp.initial <- model$getLogProb(s.nodes)
        if(n.y.mark>0){
          lp.initial <- lp.initial+model$getLogProb(y.mark.nodes)
        }
        if(n.y.mID>0){
          lp.initial <- lp.initial+model$getLogProb(y.mID.nodes)
        }
        if(n.y.mnoID>0){
          lp.initial <- lp.initial+model$getLogProb(y.mnoID.nodes)
        }
        if(n.y.um>0){
          lp.initial <- lp.initial+model$getLogProb(y.um.nodes)
        }
        if(n.y.unk>0){
          lp.initial <- lp.initial+model$getLogProb(y.unk.nodes)
        }
        
        #save current aggregate rates and subtract this individual's old contribution
        bigLam.marked.initial <- model$bigLam.marked
        bigLam.unmarked.initial <- model$bigLam.unmarked
        bigLam.marked.proposed <- bigLam.marked.initial
        bigLam.unmarked.proposed <- bigLam.unmarked.initial
        if(J.sight[g]>0&K.sight[g]>0){
          #secondary occasions when focal individual is marked
          if(n.k.marked>0){
            for(k2 in 1:n.k.marked){
              k <- k.marked[k2]
              for(j in 1:J.sight[g]){
                bigLam.old <- bigLam.marked.proposed[g,j,k]
                bigLam.marked.proposed[g,j,k] <- bigLam.old-model$lam[i,g,j]
                #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                if(bigLam.old>0&bigLam.marked.proposed[g,j,k]<=1e-12*bigLam.old){
                  bigLam.marked.proposed[g,j,k] <- 0
                  for(ii in 1:n.marked.all){
                    if(ii!=i&mark.states.all[ii,g,k]==1){
                      bigLam.marked.proposed[g,j,k] <- bigLam.marked.proposed[g,j,k]+model$lam[ii,g,j]
                    }
                  }
                }
                if(bigLam.marked.proposed[g,j,k]<0){
                  bigLam.marked.proposed[g,j,k] <- 0
                }
              }
            }
          }
          #secondary occasions when focal individual is unmarked
          if(n.k.unmarked>0){
            for(k2 in 1:n.k.unmarked){
              k <- k.unmarked[k2]3
              for(j in 1:J.sight[g]){
                bigLam.old <- bigLam.unmarked.proposed[g,j,k]
                bigLam.unmarked.proposed[g,j,k] <- bigLam.old-model$lam[i,g,j]
                #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                if(bigLam.old>0&bigLam.unmarked.proposed[g,j,k]<=1e-12*bigLam.old){
                  bigLam.unmarked.proposed[g,j,k] <- 0
                  for(ii in 1:M){
                    if(ii!=i&mark.states.all[ii,g,k]==0){
                      bigLam.unmarked.proposed[g,j,k] <- bigLam.unmarked.proposed[g,j,k]+model$lam[ii,g,j]
                    }
                  }
                }
                if(bigLam.unmarked.proposed[g,j,k]<0){
                  bigLam.unmarked.proposed[g,j,k] <- 0
                }
              }
            }
          }
        }
        
        #update proposed s and direct focal-individual observation model
        model$s[i,g,1:2] <<- s.cand
        lp.proposed <- model$calculate(s.nodes)
        #marking: s -> pd -> y.mark
        if(n.pd>0){
          model$calculate(pd.nodes)
        }
        if(n.y.mark>0){
          lp.proposed <- lp.proposed+model$calculate(y.mark.nodes)
        }
        #identified sightings: s -> lam -> y.mID; y.mID is retained only for marked k
        if(n.lam>0){
          model$calculate(lam.nodes)
        }
        if(n.y.mID>0){
          lp.proposed <- lp.proposed+model$calculate(y.mID.nodes)
        }
        
        #add this individual's new lam contribution back into the appropriate aggregate pool
        if(J.sight[g]>0&K.sight[g]>0){
          if(n.k.marked>0){
            for(k2 in 1:n.k.marked){
              k <- k.marked[k2]
              bigLam.marked.proposed[g,1:J.sight[g],k] <- bigLam.marked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
            }
          }
          if(n.k.unmarked>0){
            for(k2 in 1:n.k.unmarked){
              k <- k.unmarked[k2]
              bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
            }
          }
          model$bigLam.marked <<- bigLam.marked.proposed
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
        }
        
        #update affected aggregate sighting likelihoods
        if(n.lam.mnoID>0){
          model$calculate(lam.mnoID.nodes)
        }
        if(n.y.mnoID>0){
          lp.proposed <- lp.proposed+model$calculate(y.mnoID.nodes)
        }
        if(n.lam.um>0){
          model$calculate(lam.um.nodes)
        }
        if(n.y.um>0){
          lp.proposed <- lp.proposed+model$calculate(y.um.nodes)
        }
        if(n.lam.unk>0){
          model$calculate(lam.unk.nodes)
        }
        if(n.y.unk>0){
          lp.proposed <- lp.proposed+model$calculate(y.unk.nodes)
        }
        
        log_MH_ratio <- lp.proposed-lp.initial
        accept <- decide(log_MH_ratio)
        if(accept){
          copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
        }else{
          copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
        }
        if(adaptive){
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


#sSampler 2 updates s when z.super=1 and z=0
sSampler2 <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    # calcNodes <- model$getDependencies(target)
    calcNodes <- control$calcNodes
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
    z.super <- model$z.super[i]
    z <- model$z[i,g]
    if(z.super==1&z==0){
      s.cand <- c(rnorm(1,model$s[i,g,1],scale), rnorm(1,model$s[i,g,2],scale))
      inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        model_lp_initial <- model$getLogProb(calcNodes)
        model$s[i,g,1:2] <<- s.cand
        model_lp_proposed <- model$calculate(calcNodes)
        log_MH_ratio <- model_lp_proposed - model_lp_initial
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


#This one is designed to better jump over gaps in habitat mask
sSampler3 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    jump.multiplier <- control$jump.multiplier
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    K.sight <- control$K.sight
    n.marked.all <- control$n.marked.all
    mark.states <- control$mark.states
    mark.states.all <- control$mark.states.all
    M <- dim(mark.states.all)[1]
    s.nodes <- control$s.nodes
    calcNodes <- control$calcNodes
    
    #Explicit observation-model node sets for this primary occasion.
    #Secondary occasions are separated by the focal individual's fixed mark state.
    pd.nodes <- y.mark.nodes <- character(0)
    lam.nodes <- y.mID.nodes <- character(0)
    lam.mnoID.nodes <- y.mnoID.nodes <- character(0)
    lam.um.nodes <- y.um.nodes <- character(0)
    lam.unk.nodes <- y.unk.nodes <- character(0)
    k.marked <- k.unmarked <- integer(0)
    
    if(J.mark[g]>0){
      pd.nodes <- model$expandNodeNames(paste0("pd[",i,",",g,",1:",J.mark[g],"]"))
      y.mark.nodes <- model$expandNodeNames(paste0("y.mark[",i,",",g,",1:",J.mark[g],"]"))
    }
    if(J.sight[g]>0){
      lam.nodes <- model$expandNodeNames(paste0("lam[",i,",",g,",1:",J.sight[g],"]"))
      if(K.sight[g]>0){
        for(k.setup in 1:K.sight[g]){
          if(mark.states[g,k.setup]==1){
            k.marked <- c(k.marked,k.setup)
            if(i<=n.marked.all){
              y.mID.nodes <- c(y.mID.nodes,
                               model$expandNodeNames(paste0("y.mID[",i,",",g,",1:",J.sight[g],",",k.setup,"]")))
            }
            lam.mnoID.nodes <- c(lam.mnoID.nodes,
                                 model$expandNodeNames(paste0("lam.mnoID[",g,",1:",J.sight[g],",",k.setup,"]")))
            y.mnoID.nodes <- c(y.mnoID.nodes,
                               model$expandNodeNames(paste0("y.mnoID[",g,",1:",J.sight[g],",",k.setup,"]")))
          }else{
            k.unmarked <- c(k.unmarked,k.setup)
            lam.um.nodes <- c(lam.um.nodes,
                              model$expandNodeNames(paste0("lam.um[",g,",1:",J.sight[g],",",k.setup,"]")))
            y.um.nodes <- c(y.um.nodes,
                            model$expandNodeNames(paste0("y.um[",g,",1:",J.sight[g],",",k.setup,"]")))
          }
          #unknown-mark-status sightings always depend on the focal individual through
          #either bigLam.marked or bigLam.unmarked
          lam.unk.nodes <- c(lam.unk.nodes,
                             model$expandNodeNames(paste0("lam.unk[",g,",1:",J.sight[g],",",k.setup,"]")))
          y.unk.nodes <- c(y.unk.nodes,
                           model$expandNodeNames(paste0("y.unk[",g,",1:",J.sight[g],",",k.setup,"]")))
        }
      }
    }
    n.pd <- length(pd.nodes)
    n.y.mark <- length(y.mark.nodes)
    n.lam <- length(lam.nodes)
    n.y.mID <- length(y.mID.nodes)
    n.lam.mnoID <- length(lam.mnoID.nodes)
    n.y.mnoID <- length(y.mnoID.nodes)
    n.lam.um <- length(lam.um.nodes)
    n.y.um <- length(y.um.nodes)
    n.lam.unk <- length(lam.unk.nodes)
    n.y.unk <- length(y.unk.nodes)
    n.k.marked <- length(k.marked)
    n.k.unmarked <- length(k.unmarked)
  },
  run = function(){
    if(model$z.super[i]==1){
      if(g==1){
        scale.jump <- jump.multiplier * model$sigma.move.int[1]
      }else{
        scale.jump <- jump.multiplier * model$sigma.move.int[g-1]
      }
      s.cand <- c(rnorm(1,model$s[i,g,1],scale.jump),
                  rnorm(1,model$s[i,g,2],scale.jump))
      inbox <- s.cand[1] < xlim[2] & s.cand[1] > xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        if(model$z[i,g]==0){
          model_lp_initial <- model$getLogProb(calcNodes)
          model$s[i,g,1:2] <<- s.cand
          model_lp_proposed <- model$calculate(calcNodes)
          log_MH_ratio <- model_lp_proposed - model_lp_initial
          accept <- decide(log_MH_ratio)
          if(accept){
            copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
          }else{
            copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
          }
        }else{
          #initial log probability for movement/telemetry terms and affected observation likelihoods
          lp.initial <- model$getLogProb(s.nodes)
          if(n.y.mark>0){
            lp.initial <- lp.initial+model$getLogProb(y.mark.nodes)
          }
          if(n.y.mID>0){
            lp.initial <- lp.initial+model$getLogProb(y.mID.nodes)
          }
          if(n.y.mnoID>0){
            lp.initial <- lp.initial+model$getLogProb(y.mnoID.nodes)
          }
          if(n.y.um>0){
            lp.initial <- lp.initial+model$getLogProb(y.um.nodes)
          }
          if(n.y.unk>0){
            lp.initial <- lp.initial+model$getLogProb(y.unk.nodes)
          }
          
          #save current aggregate rates and subtract this individual's old contribution
          bigLam.marked.initial <- model$bigLam.marked
          bigLam.unmarked.initial <- model$bigLam.unmarked
          bigLam.marked.proposed <- bigLam.marked.initial
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          if(J.sight[g]>0&K.sight[g]>0){
            #secondary occasions when focal individual is marked
            if(n.k.marked>0){
              for(k2 in 1:n.k.marked){
                k <- k.marked[k2]
                for(j in 1:J.sight[g]){
                  bigLam.old <- bigLam.marked.proposed[g,j,k]
                  bigLam.marked.proposed[g,j,k] <- bigLam.old-model$lam[i,g,j]
                  #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                  if(bigLam.old>0&bigLam.marked.proposed[g,j,k]<=1e-12*bigLam.old){
                    bigLam.marked.proposed[g,j,k] <- 0
                    for(ii in 1:n.marked.all){
                      if(ii!=i&mark.states.all[ii,g,k]==1){
                        bigLam.marked.proposed[g,j,k] <- bigLam.marked.proposed[g,j,k]+model$lam[ii,g,j]
                      }
                    }
                  }
                  if(bigLam.marked.proposed[g,j,k]<0){
                    bigLam.marked.proposed[g,j,k] <- 0
                  }
                }
              }
            }
            #secondary occasions when focal individual is unmarked
            if(n.k.unmarked>0){
              for(k2 in 1:n.k.unmarked){
                k <- k.unmarked[k2]
                for(j in 1:J.sight[g]){
                  bigLam.old <- bigLam.unmarked.proposed[g,j,k]
                  bigLam.unmarked.proposed[g,j,k] <- bigLam.old-model$lam[i,g,j]
                  #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                  if(bigLam.old>0&bigLam.unmarked.proposed[g,j,k]<=1e-12*bigLam.old){
                    bigLam.unmarked.proposed[g,j,k] <- 0
                    for(ii in 1:M){
                      if(ii!=i&mark.states.all[ii,g,k]==0){
                        bigLam.unmarked.proposed[g,j,k] <- bigLam.unmarked.proposed[g,j,k]+model$lam[ii,g,j]
                      }
                    }
                  }
                  if(bigLam.unmarked.proposed[g,j,k]<0){
                    bigLam.unmarked.proposed[g,j,k] <- 0
                  }
                }
              }
            }
          }
          
          #update proposed s and direct focal-individual observation model
          model$s[i,g,1:2] <<- s.cand
          lp.proposed <- model$calculate(s.nodes)
          #marking: s -> pd -> y.mark
          if(n.pd>0){
            model$calculate(pd.nodes)
          }
          if(n.y.mark>0){
            lp.proposed <- lp.proposed+model$calculate(y.mark.nodes)
          }
          #identified sightings: s -> lam -> y.mID; y.mID is retained only for marked k
          if(n.lam>0){
            model$calculate(lam.nodes)
          }
          if(n.y.mID>0){
            lp.proposed <- lp.proposed+model$calculate(y.mID.nodes)
          }
          
          #add this individual's new lam contribution back into the appropriate aggregate pool
          if(J.sight[g]>0&K.sight[g]>0){
            if(n.k.marked>0){
              for(k2 in 1:n.k.marked){
                k <- k.marked[k2]
                bigLam.marked.proposed[g,1:J.sight[g],k] <- bigLam.marked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
              }
            }
            if(n.k.unmarked>0){
              for(k2 in 1:n.k.unmarked){
                k <- k.unmarked[k2]
                bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
              }
            }
            model$bigLam.marked <<- bigLam.marked.proposed
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
          }
          
          #update affected aggregate sighting likelihoods
          if(n.lam.mnoID>0){
            model$calculate(lam.mnoID.nodes)
          }
          if(n.y.mnoID>0){
            lp.proposed <- lp.proposed+model$calculate(y.mnoID.nodes)
          }
          if(n.lam.um>0){
            model$calculate(lam.um.nodes)
          }
          if(n.y.um>0){
            lp.proposed <- lp.proposed+model$calculate(y.um.nodes)
          }
          if(n.lam.unk>0){
            model$calculate(lam.unk.nodes)
          }
          if(n.y.unk>0){
            lp.proposed <- lp.proposed+model$calculate(y.unk.nodes)
          }
          
          log_MH_ratio <- lp.proposed-lp.initial
          accept <- decide(log_MH_ratio)
          if(accept){
            copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
          }else{
            copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
          }
        }
      }
    }
  },
  methods = list(
    reset = function() {}
  )
)

