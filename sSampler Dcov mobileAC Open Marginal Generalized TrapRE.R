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
    n.marked.all <- control$n.marked.all
    mark.states <- control$mark.states
    s.nodes <- control$s.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    y.mark.nodes <- control$y.mark.nodes
    y.mID.nodes <- control$y.mID.nodes
    y.mnoID.nodes <- control$y.mnoID.nodes
    y.um.nodes <- control$y.um.nodes
    y.unk.nodes <- control$y.unk.nodes
    trap.RE.nodes <- control$trap.RE.nodes
    lam.mnoID.nodes <- control$lam.mnoID.nodes
    lam.um.nodes <- control$lam.um.nodes
    lam.unk.nodes <- control$lam.unk.nodes
    calcNodes <- control$calcNodes
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
  run <- function(){
    z.super <- model$z.super[i]
    z <- model$z[i,g]
    if(z.super==1&z==1){
      s.cand <- c(rnorm(1,model$s[i,g,1],scale),rnorm(1,model$s[i,g,2],scale))
      inbox <- s.cand[1]<xlim[2] & s.cand[1]>xlim[1] & s.cand[2]<ylim[2] & s.cand[2]>ylim[1]
      if(inbox){
        lp.initial.s <- model$getLogProb(s.nodes)
        lp.initial.y.mark <- model$getLogProb(y.mark.nodes)
        lp.initial.y.um <- model$getLogProb(y.um.nodes)
        lp.initial.y.unk <- model$getLogProb(y.unk.nodes)
        lp.initial.trap.RE <- model$getLogProb(trap.RE.nodes)
        if(i<=n.marked.all){
          lp.initial.y.mID <- model$getLogProb(y.mID.nodes)
          lp.initial.y.mnoID <- model$getLogProb(y.mnoID.nodes)
          bigLam.marked.initial <- model$bigLam.marked
          bigLam.unmarked.initial <- model$bigLam.unmarked
          model$s[i,g,1:2] <<- s.cand
          lp.proposed.s <- model$calculate(s.nodes)
          bigLam.marked.proposed <- bigLam.marked.initial
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          if(J.sight[g]>0){
            bigLam.marked.proposed[g,1:J.sight[g]] <- bigLam.marked.proposed[g,1:J.sight[g]] - model$lam[i,g,1:J.sight[g]]*mark.states[g]
            bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] - model$lam[i,g,1:J.sight[g]]*(1-mark.states[g])
            for(j in 1:J.sight[g]){
              if(bigLam.marked.proposed[g,j]<0){
                bigLam.marked.proposed[g,j] <- 0
              }
              if(bigLam.unmarked.proposed[g,j]<0){
                bigLam.unmarked.proposed[g,j] <- 0
              }
            }
          }
          if(J.mark[g]>0){
            model$calculate(pd.nodes)
          }
          if(J.sight[g]>0){
            model$calculate(lam.nodes)
            bigLam.marked.proposed[g,1:J.sight[g]] <- bigLam.marked.proposed[g,1:J.sight[g]] + model$lam[i,g,1:J.sight[g]]*mark.states[g]
            bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] + model$lam[i,g,1:J.sight[g]]*(1-mark.states[g])
            model$bigLam.marked <<- bigLam.marked.proposed
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            model$calculate(lam.mnoID.nodes)
            model$calculate(lam.um.nodes)
            model$calculate(lam.unk.nodes)
          }
          lp.proposed.y.mark <- model$calculate(y.mark.nodes)
          lp.proposed.y.mID <- model$calculate(y.mID.nodes)
          lp.proposed.y.mnoID <- model$calculate(y.mnoID.nodes)
          lp.proposed.y.um <- model$calculate(y.um.nodes)
          lp.proposed.y.unk <- model$calculate(y.unk.nodes)
          lp.proposed.trap.RE <- model$calculate(trap.RE.nodes)
          lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.mID + lp.initial.y.mnoID +
            lp.initial.y.um + lp.initial.y.unk + lp.initial.trap.RE
          lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.mID + lp.proposed.y.mnoID +
            lp.proposed.y.um + lp.proposed.y.unk + lp.proposed.trap.RE
        }else{
          bigLam.unmarked.initial <- model$bigLam.unmarked
          model$s[i,g,1:2] <<- s.cand
          lp.proposed.s <- model$calculate(s.nodes)
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          if(J.sight[g]>0){
            bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] - model$lam[i,g,1:J.sight[g]]
            for(j in 1:J.sight[g]){
              if(bigLam.unmarked.proposed[g,j]<0){
                bigLam.unmarked.proposed[g,j] <- 0
              }
            }
          }
          if(J.mark[g]>0){
            model$calculate(pd.nodes)
          }
          if(J.sight[g]>0){
            model$calculate(lam.nodes)
            bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] + model$lam[i,g,1:J.sight[g]]
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            model$calculate(lam.um.nodes)
            model$calculate(lam.unk.nodes)
          }
          lp.proposed.y.mark <- model$calculate(y.mark.nodes)
          lp.proposed.y.um <- model$calculate(y.um.nodes)
          lp.proposed.y.unk <- model$calculate(y.unk.nodes)
          lp.proposed.trap.RE <- model$calculate(trap.RE.nodes)
          lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk + lp.initial.trap.RE
          lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk + lp.proposed.trap.RE
        }
        log_MH_ratio <- lp.proposed - lp.initial
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
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    jump.multiplier <- control$jump.multiplier
    sig.move.fixed <- control$sig.move.fixed
    calcNodes <- control$calcNodes
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    n.marked.all <- control$n.marked.all
    mark.states <- control$mark.states
    s.nodes <- control$s.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    y.mark.nodes <- control$y.mark.nodes
    y.mID.nodes <- control$y.mID.nodes
    y.mnoID.nodes <- control$y.mnoID.nodes
    y.um.nodes <- control$y.um.nodes
    y.unk.nodes <- control$y.unk.nodes
    lam.mnoID.nodes <- control$lam.mnoID.nodes
    lam.um.nodes <- control$lam.um.nodes
    lam.unk.nodes <- control$lam.unk.nodes
    trap.RE.nodes <- control$trap.RE.nodes
  },
  run = function(){
    if(model$z.super[i]==1){
      if(sig.move.fixed==TRUE){
        scale.jump <- jump.multiplier*model$sigma.move[1]
      }else{
        scale.jump <- jump.multiplier*model$sigma.move[i]
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
          lp.initial.s <- model$getLogProb(s.nodes)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes)
          lp.initial.y.um <- model$getLogProb(y.um.nodes)
          lp.initial.y.unk <- model$getLogProb(y.unk.nodes)
          lp.initial.trap.RE <- model$getLogProb(trap.RE.nodes)
          if(i<=n.marked.all){
            lp.initial.y.mID <- model$getLogProb(y.mID.nodes)
            lp.initial.y.mnoID <- model$getLogProb(y.mnoID.nodes)
            bigLam.marked.initial <- model$bigLam.marked
            bigLam.unmarked.initial <- model$bigLam.unmarked
            model$s[i,g,1:2] <<- s.cand
            lp.proposed.s <- model$calculate(s.nodes)
            bigLam.marked.proposed <- bigLam.marked.initial
            bigLam.unmarked.proposed <- bigLam.unmarked.initial
            if(J.sight[g]>0){
              bigLam.marked.proposed[g,1:J.sight[g]] <- bigLam.marked.proposed[g,1:J.sight[g]] - model$lam[i,g,1:J.sight[g]]*mark.states[g]
              bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] - model$lam[i,g,1:J.sight[g]]*(1-mark.states[g])
              for(j in 1:J.sight[g]){
                if(bigLam.marked.proposed[g,j]<0){
                  bigLam.marked.proposed[g,j] <- 0
                }
                if(bigLam.unmarked.proposed[g,j]<0){
                  bigLam.unmarked.proposed[g,j] <- 0
                }
              }
            }
            if(J.mark[g]>0){
              model$calculate(pd.nodes)
            }
            if(J.sight[g]>0){
              model$calculate(lam.nodes)
              bigLam.marked.proposed[g,1:J.sight[g]] <- bigLam.marked.proposed[g,1:J.sight[g]] + model$lam[i,g,1:J.sight[g]]*mark.states[g]
              bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] + model$lam[i,g,1:J.sight[g]]*(1-mark.states[g])
              model$bigLam.marked <<- bigLam.marked.proposed
              model$bigLam.unmarked <<- bigLam.unmarked.proposed
              model$calculate(lam.mnoID.nodes)
              model$calculate(lam.um.nodes)
              model$calculate(lam.unk.nodes)
            }
            lp.proposed.y.mark <- model$calculate(y.mark.nodes)
            lp.proposed.y.mID <- model$calculate(y.mID.nodes)
            lp.proposed.y.mnoID <- model$calculate(y.mnoID.nodes)
            lp.proposed.y.um <- model$calculate(y.um.nodes)
            lp.proposed.y.unk <- model$calculate(y.unk.nodes)
            lp.proposed.trap.RE <- model$calculate(trap.RE.nodes)
            lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.mID + lp.initial.y.mnoID +
              lp.initial.y.um + lp.initial.y.unk + lp.initial.trap.RE
            lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.mID + lp.proposed.y.mnoID +
              lp.proposed.y.um + lp.proposed.y.unk + lp.proposed.trap.RE
          }else{
            bigLam.unmarked.initial <- model$bigLam.unmarked
            model$s[i,g,1:2] <<- s.cand
            lp.proposed.s <- model$calculate(s.nodes)
            bigLam.unmarked.proposed <- bigLam.unmarked.initial
            if(J.sight[g]>0){
              bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] - model$lam[i,g,1:J.sight[g]]
              for(j in 1:J.sight[g]){
                if(bigLam.unmarked.proposed[g,j]<0){
                  bigLam.unmarked.proposed[g,j] <- 0
                }
              }
            }
            if(J.mark[g]>0){
              model$calculate(pd.nodes)
            }
            if(J.sight[g]>0){
              model$calculate(lam.nodes)
              bigLam.unmarked.proposed[g,1:J.sight[g]] <- bigLam.unmarked.proposed[g,1:J.sight[g]] + model$lam[i,g,1:J.sight[g]]
              model$bigLam.unmarked <<- bigLam.unmarked.proposed
              model$calculate(lam.um.nodes)
              model$calculate(lam.unk.nodes)
            }
            lp.proposed.y.mark <- model$calculate(y.mark.nodes)
            lp.proposed.y.um <- model$calculate(y.um.nodes)
            lp.proposed.y.unk <- model$calculate(y.unk.nodes)
            lp.proposed.trap.RE <- model$calculate(trap.RE.nodes)
            lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk + lp.initial.trap.RE
            lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk + lp.proposed.trap.RE
          }
          log_MH_ratio <- lp.proposed - lp.initial
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