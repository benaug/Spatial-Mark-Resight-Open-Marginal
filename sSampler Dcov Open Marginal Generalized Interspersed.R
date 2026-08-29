sSamplerDcov <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    i <- control$i
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    K.sight <- control$K.sight
    sight.g <- control$sight.g
    n.sight.g <- control$n.sight.g
    res <- control$res
    xlim <- control$xlim
    ylim <- control$ylim
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    n.marked.all <- control$n.marked.all
    n.primary <- control$n.primary
    mark.states <- control$mark.states
    mark.states.all <- control$mark.states.all
    M <- dim(mark.states.all)[1]
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    #full dependency set retained only for saving/restoring model state after accept/reject
    calcNodes <- model$getDependencies(target)
    loc.nodes <- control$loc.nodes
    s.nodes <- c(model$expandNodeNames(paste("s[",i,",",1:2,"]")),
                 model$expandNodeNames(paste("s.cell[",i,"]")),
                 model$expandNodeNames(paste("dummy.data[",i,"]")))
    #if we have telemetry for this individual, add locs to s.nodes
    if(length(loc.nodes)>0){
      s.nodes <- c(s.nodes,loc.nodes)
    }
    
    #Explicit observation-model node sets. Nodes are grouped by primary occasion
    #with start/count vectors so run() can skip occasions where z[i,g]==0.
    pd.nodes <- y.mark.nodes <- character(0)
    lam.nodes <- y.mID.nodes <- character(0)
    lam.mnoID.nodes <- y.mnoID.nodes <- character(0)
    lam.um.nodes <- y.um.nodes <- character(0)
    lam.unk.nodes <- y.unk.nodes <- character(0)
    
    pd.node.start <- integer(n.primary)
    pd.node.n <- integer(n.primary)
    y.mark.node.start <- integer(n.primary)
    y.mark.node.n <- integer(n.primary)
    lam.node.start <- integer(n.primary)
    lam.node.n <- integer(n.primary)
    y.mID.node.start <- integer(n.primary)
    y.mID.node.n <- integer(n.primary)
    lam.mnoID.node.start <- integer(n.primary)
    lam.mnoID.node.n <- integer(n.primary)
    y.mnoID.node.start <- integer(n.primary)
    y.mnoID.node.n <- integer(n.primary)
    lam.um.node.start <- integer(n.primary)
    lam.um.node.n <- integer(n.primary)
    y.um.node.start <- integer(n.primary)
    y.um.node.n <- integer(n.primary)
    lam.unk.node.start <- integer(n.primary)
    lam.unk.node.n <- integer(n.primary)
    y.unk.node.start <- integer(n.primary)
    y.unk.node.n <- integer(n.primary)
    
    #Secondary-occasion indices by focal mark state. This removes mark.states[g,k]
    #tests from run() and allows arbitrary within-primary mark loss and remarking.
    k.marked <- k.unmarked <- integer(0)
    k.marked.start <- integer(n.primary)
    k.marked.n <- integer(n.primary)
    k.unmarked.start <- integer(n.primary)
    k.unmarked.n <- integer(n.primary)
    for(g.setup in 1:n.primary){
      pd.node.start[g.setup] <- length(pd.nodes)+1
      y.mark.node.start[g.setup] <- length(y.mark.nodes)+1
      lam.node.start[g.setup] <- length(lam.nodes)+1
      y.mID.node.start[g.setup] <- length(y.mID.nodes)+1
      lam.mnoID.node.start[g.setup] <- length(lam.mnoID.nodes)+1
      y.mnoID.node.start[g.setup] <- length(y.mnoID.nodes)+1
      lam.um.node.start[g.setup] <- length(lam.um.nodes)+1
      y.um.node.start[g.setup] <- length(y.um.nodes)+1
      lam.unk.node.start[g.setup] <- length(lam.unk.nodes)+1
      y.unk.node.start[g.setup] <- length(y.unk.nodes)+1
      k.marked.start[g.setup] <- length(k.marked)+1
      k.unmarked.start[g.setup] <- length(k.unmarked)+1
      if(J.mark[g.setup]>0){
        pd.nodes <- c(pd.nodes,
                      model$expandNodeNames(paste0("pd[",i,",",g.setup,",1:",J.mark[g.setup],"]")))
        y.mark.nodes <- c(y.mark.nodes,
                          model$expandNodeNames(paste0("y.mark[",i,",",g.setup,",1:",J.mark[g.setup],"]")))
      }
      if(J.sight[g.setup]>0&K.sight[g.setup]>0){
        lam.nodes <- c(lam.nodes,
                       model$expandNodeNames(paste0("lam[",i,",",g.setup,",1:",J.sight[g.setup],"]")))
        for(k.setup in 1:K.sight[g.setup]){
          if(mark.states[g.setup,k.setup]==1){
            k.marked <- c(k.marked,k.setup)
            if(i<=n.marked.all){
              y.mID.nodes <- c(y.mID.nodes,
                               model$expandNodeNames(paste0("y.mID[",i,",",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
            }
            lam.mnoID.nodes <- c(lam.mnoID.nodes,
                                 model$expandNodeNames(paste0("lam.mnoID[",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
            y.mnoID.nodes <- c(y.mnoID.nodes,
                               model$expandNodeNames(paste0("y.mnoID[",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
          }else{
            k.unmarked <- c(k.unmarked,k.setup)
            lam.um.nodes <- c(lam.um.nodes,
                              model$expandNodeNames(paste0("lam.um[",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
            y.um.nodes <- c(y.um.nodes,
                            model$expandNodeNames(paste0("y.um[",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
          }
          #unknown-mark-status sightings always depend on the focal individual through
          #either bigLam.marked or bigLam.unmarked
          lam.unk.nodes <- c(lam.unk.nodes,
                             model$expandNodeNames(paste0("lam.unk[",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
          y.unk.nodes <- c(y.unk.nodes,
                           model$expandNodeNames(paste0("y.unk[",g.setup,",1:",J.sight[g.setup],",",k.setup,"]")))
        }
      }
      pd.node.n[g.setup] <- length(pd.nodes)-pd.node.start[g.setup]+1
      y.mark.node.n[g.setup] <- length(y.mark.nodes)-y.mark.node.start[g.setup]+1
      lam.node.n[g.setup] <- length(lam.nodes)-lam.node.start[g.setup]+1
      y.mID.node.n[g.setup] <- length(y.mID.nodes)-y.mID.node.start[g.setup]+1
      lam.mnoID.node.n[g.setup] <- length(lam.mnoID.nodes)-lam.mnoID.node.start[g.setup]+1
      y.mnoID.node.n[g.setup] <- length(y.mnoID.nodes)-y.mnoID.node.start[g.setup]+1
      lam.um.node.n[g.setup] <- length(lam.um.nodes)-lam.um.node.start[g.setup]+1
      y.um.node.n[g.setup] <- length(y.um.nodes)-y.um.node.start[g.setup]+1
      lam.unk.node.n[g.setup] <- length(lam.unk.nodes)-lam.unk.node.start[g.setup]+1
      y.unk.node.n[g.setup] <- length(y.unk.nodes)-y.unk.node.start[g.setup]+1
      k.marked.n[g.setup] <- length(k.marked)-k.marked.start[g.setup]+1
      k.unmarked.n[g.setup] <- length(k.unmarked)-k.unmarked.start[g.setup]+1
    }
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
  run = function(){
    z.super <- model$z.super[i]
    if(z.super==0){ #propose from prior
      #propose new cell
      model$s.cell[i] <<- rcat(1,model$pi.cell)
      #propose x and y in new cell
      s.cell.x <- model$s.cell[i]%%n.cells.x
      s.cell.y <- floor(model$s.cell[i]/n.cells.x)+1
      if(s.cell.x==0){
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y-1
      }
      xlim.cell <- c(s.cell.x-1,s.cell.x)*res
      ylim.cell <- c(s.cell.y-1,s.cell.y)*res
      model$s[i,1:2] <<- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
      model$calculate(s.nodes)
      copy(from = model, to = mvSaved, row = 1, nodes = s.nodes, logProb = TRUE)
    }else{ #MH
      s.cand <- c(rnorm(1,model$s[i,1],scale), rnorm(1,model$s[i,2],scale))
      inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        #initial log probability for s/density/telemetry terms and observation likelihoods
        #only in primary occasions where this individual is alive
        lp.initial <- model$getLogProb(s.nodes)
        for(g in 1:n.primary){
          if(model$z[i,g]==1){
            if(y.mark.node.n[g]>0){
              node.start <- y.mark.node.start[g]
              node.end <- node.start+y.mark.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(y.mark.nodes[node.start:node.end])
            }
            if(y.mID.node.n[g]>0){
              node.start <- y.mID.node.start[g]
              node.end <- node.start+y.mID.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(y.mID.nodes[node.start:node.end])
            }
            if(y.mnoID.node.n[g]>0){
              node.start <- y.mnoID.node.start[g]
              node.end <- node.start+y.mnoID.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(y.mnoID.nodes[node.start:node.end])
            }
            if(y.um.node.n[g]>0){
              node.start <- y.um.node.start[g]
              node.end <- node.start+y.um.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(y.um.nodes[node.start:node.end])
            }
            if(y.unk.node.n[g]>0){
              node.start <- y.unk.node.start[g]
              node.end <- node.start+y.unk.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(y.unk.nodes[node.start:node.end])
            }
          }
        }
        
        #save current aggregate rates and subtract this individual's old contribution
        bigLam.marked.initial <- model$bigLam.marked
        bigLam.unmarked.initial <- model$bigLam.unmarked
        bigLam.marked.proposed <- bigLam.marked.initial
        bigLam.unmarked.proposed <- bigLam.unmarked.initial
        for(g2 in 1:n.sight.g){
          g <- sight.g[g2]
          if(model$z[i,g]==1){ #z.super always 1 here
            #secondary occasions when focal individual is marked
            if(k.marked.n[g]>0){
              k.start <- k.marked.start[g]
              k.end <- k.start+k.marked.n[g]-1
              for(k2 in k.start:k.end){
                k <- k.marked[k2]
                for(j in 1:J.sight[g]){
                  bigLam.old <- bigLam.marked.proposed[g,j,k]
                  bigLam.marked.proposed[g,j,k] <- bigLam.old-model$lam[i,g,j]
                  #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                  if(bigLam.old>0&bigLam.marked.proposed[g,j,k]<=1e-12*bigLam.old){
                    bigLam.marked.proposed[g,j,k] <- 0
                    for(ii in 1:M){
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
            if(k.unmarked.n[g]>0){
              k.start <- k.unmarked.start[g]
              k.end <- k.start+k.unmarked.n[g]-1
              for(k2 in k.start:k.end){
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
        }
        
        #update proposed s, then the direct focal-individual observation model
        model$s[i, 1:2] <<- s.cand
        lp.proposed <- model$calculate(s.nodes)
        for(g in 1:n.primary){
          if(model$z[i,g]==1){
            #marking: s -> pd -> y.mark
            if(pd.node.n[g]>0){
              node.start <- pd.node.start[g]
              node.end <- node.start+pd.node.n[g]-1
              model$calculate(pd.nodes[node.start:node.end])
            }
            if(y.mark.node.n[g]>0){
              node.start <- y.mark.node.start[g]
              node.end <- node.start+y.mark.node.n[g]-1
              lp.proposed <- lp.proposed+model$calculate(y.mark.nodes[node.start:node.end])
            }
            #identified sightings: s -> lam -> y.mID; y.mID nodes are retained only for marked k
            if(lam.node.n[g]>0){
              node.start <- lam.node.start[g]
              node.end <- node.start+lam.node.n[g]-1
              model$calculate(lam.nodes[node.start:node.end])
            }
            if(y.mID.node.n[g]>0){
              node.start <- y.mID.node.start[g]
              node.end <- node.start+y.mID.node.n[g]-1
              lp.proposed <- lp.proposed+model$calculate(y.mID.nodes[node.start:node.end])
            }
          }
        }
        
        #add this individual's new lam contribution back into the appropriate aggregate pool
        for(g2 in 1:n.sight.g){
          g <- sight.g[g2]
          if(model$z[i,g]==1){
            if(k.marked.n[g]>0){
              k.start <- k.marked.start[g]
              k.end <- k.start+k.marked.n[g]-1
              for(k2 in k.start:k.end){
                k <- k.marked[k2]
                bigLam.marked.proposed[g,1:J.sight[g],k] <- bigLam.marked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
              }
            }
            if(k.unmarked.n[g]>0){
              k.start <- k.unmarked.start[g]
              k.end <- k.start+k.unmarked.n[g]-1
              for(k2 in k.start:k.end){
                k <- k.unmarked[k2]
                bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
              }
            }
          }
        }
        model$bigLam.marked <<- bigLam.marked.proposed
        model$bigLam.unmarked <<- bigLam.unmarked.proposed
        
        #update aggregate sighting likelihoods whose rates changed for this individual
        for(g in 1:n.primary){
          if(model$z[i,g]==1){
            #marked no-ID sightings only change on secondary occasions when focal individual is marked
            if(lam.mnoID.node.n[g]>0){
              node.start <- lam.mnoID.node.start[g]
              node.end <- node.start+lam.mnoID.node.n[g]-1
              model$calculate(lam.mnoID.nodes[node.start:node.end])
            }
            if(y.mnoID.node.n[g]>0){
              node.start <- y.mnoID.node.start[g]
              node.end <- node.start+y.mnoID.node.n[g]-1
              lp.proposed <- lp.proposed+model$calculate(y.mnoID.nodes[node.start:node.end])
            }
            #unmarked sightings only change on secondary occasions when focal individual is unmarked
            if(lam.um.node.n[g]>0){
              node.start <- lam.um.node.start[g]
              node.end <- node.start+lam.um.node.n[g]-1
              model$calculate(lam.um.nodes[node.start:node.end])
            }
            if(y.um.node.n[g]>0){
              node.start <- y.um.node.start[g]
              node.end <- node.start+y.um.node.n[g]-1
              lp.proposed <- lp.proposed+model$calculate(y.um.nodes[node.start:node.end])
            }
            #unknown-mark-status sightings change on every secondary occasion
            if(lam.unk.node.n[g]>0){
              node.start <- lam.unk.node.start[g]
              node.end <- node.start+lam.unk.node.n[g]-1
              model$calculate(lam.unk.nodes[node.start:node.end])
            }
            if(y.unk.node.n[g]>0){
              node.start <- y.unk.node.start[g]
              node.end <- node.start+y.unk.node.n[g]-1
              lp.proposed <- lp.proposed+model$calculate(y.unk.nodes[node.start:node.end])
            }
          }
        }
        log.MH.ratio <- lp.proposed - lp.initial
        accept <- decide(log.MH.ratio)
        if(accept) {
          copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
          copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive){ #tune RW proposals when z.super=1
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