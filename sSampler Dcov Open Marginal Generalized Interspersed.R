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
    calcNodes <- model$getDependencies(target)
    loc.nodes <- control$loc.nodes
    s.nodes <- c(model$expandNodeNames(paste("s[",i,",",1:2,"]")),
                 model$expandNodeNames(paste("s.cell[",i,"]")),
                 model$expandNodeNames(paste("dummy.data[",i,"]")))
    #if we have telemetry for this individual, add locs to s.nodes
    if(length(loc.nodes)>0){
      s.nodes <- c(s.nodes,loc.nodes)
    }
    #OLD observation-specific node definitions - replaced by graph-derived primary-occasion blocks
    # #expand requested blocks to the nodes that actually exist in the model
    # pd.nodes <- model$expandNodeNames(paste("pd[",i,",1:",n.primary,",1:",max(J.mark),"]"))
    # y.mark.nodes <- model$expandNodeNames(paste("y.mark[",i,",1:",n.primary,",1:",max(J.mark),"]"))
    # y.mID.nodes <- model$expandNodeNames(paste("y.mID[",i,",1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    # y.mnoID.nodes <- model$expandNodeNames(paste("y.mnoID[1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    # y.um.nodes <- model$expandNodeNames(paste("y.um[1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    # y.unk.nodes <- model$expandNodeNames(paste("y.unk[1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    # lam.nodes <- model$expandNodeNames(paste("lam[",i,",1:",n.primary,",1:",max(J.sight),"]"))
    # lam.mnoID.nodes <- model$expandNodeNames(paste("lam.mnoID[1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    # lam.um.nodes <- model$expandNodeNames(paste("lam.um[1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    # lam.unk.nodes <- model$expandNodeNames(paste("lam.unk[1:",n.primary,",1:",max(J.sight),",1:",max(K.sight),"]"))
    
    #Group dependencies of s[i,] by primary occasion using the model graph.
    #focal.nodes contain individual observation terms (e.g., pd/y.mark and lam/y.mID).
    #aggregate.nodes contain terms downstream of bigLam (e.g., lam.mnoID/um/unk and their likelihoods).
    #bigLam itself is excluded because it is updated manually below.
    focal.nodes <- aggregate.nodes <- c()
    focal.node.start <- integer(n.primary)
    focal.node.n <- integer(n.primary)
    aggregate.node.start <- integer(n.primary)
    aggregate.node.n <- integer(n.primary)
    for(g.setup in 1:n.primary){
      z.dep <- model$getDependencies(paste0("z[",i,",",g.setup,"]"))
      local.g <- calcNodes[calcNodes %in% z.dep]
      
      #identify bigLam nodes for this primary occasion; these are updated manually
      bigLam.g <- local.g[grepl("^bigLam\\.(marked|unmarked)\\[",local.g)]
      
      #anything downstream of bigLam is an aggregate observation term
      aggregate.g <- character(0)
      if(length(bigLam.g)>0){
        aggregate.dep <- model$getDependencies(bigLam.g,self=FALSE)
        aggregate.g <- local.g[local.g %in% aggregate.dep]
      }
      
      #remaining occasion-specific dependencies are focal-individual terms
      focal.g <- local.g[!local.g %in% c(bigLam.g,aggregate.g)]
      
      focal.node.start[g.setup] <- length(focal.nodes)+1
      focal.node.n[g.setup] <- length(focal.g)
      if(length(focal.g)>0){
        focal.nodes <- c(focal.nodes,focal.g)
      }
      
      aggregate.node.start[g.setup] <- length(aggregate.nodes)+1
      aggregate.node.n[g.setup] <- length(aggregate.g)
      if(length(aggregate.g)>0){
        aggregate.nodes <- c(aggregate.nodes,aggregate.g)
      }
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
        # OLD marked/unmarked-specific full observation update - replaced by unified alive-year graph update
        # #get initial logprobs - not optimizing by considering if this is marked or unmarked i
        # lp.initial.s <- model$getLogProb(s.nodes)
        # lp.initial.y.mark <- model$getLogProb(y.mark.nodes)
        # lp.initial.y.unk <- model$getLogProb(y.unk.nodes)
        # lp.initial.y.um <- model$getLogProb(y.um.nodes)
        # if(i<=n.marked.all){ #if marked in at least one year
        #   lp.initial.y.mID <- model$getLogProb(y.mID.nodes)
        #   lp.initial.y.mnoID <- model$getLogProb(y.mnoID.nodes)
        #   #pull these out of model object
        #   bigLam.marked.initial <- model$bigLam.marked
        #   bigLam.unmarked.initial <- model$bigLam.unmarked
        #   #update proposed s
        #   model$s[i, 1:2] <<- s.cand
        #   lp.proposed.s <- model$calculate(s.nodes) #proposed logprob for s.nodes
        #   #subtract these out before calculating lam
        #   bigLam.marked.proposed <- bigLam.marked.initial
        #   bigLam.unmarked.proposed <- bigLam.unmarked.initial
        #   for(g2 in 1:n.sight.g){
        #     g <- sight.g[g2]
        #     if(model$z[i,g]==1){ #z.super always 1 here
        #       for(k in 1:K.sight[g]){
        #         bigLam.marked.proposed[g,1:J.sight[g],k] <- bigLam.marked.proposed[g,1:J.sight[g],k] - model$lam[i,g,1:J.sight[g]]*mark.states[g,k]
        #         bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] - model$lam[i,g,1:J.sight[g]]*(1-mark.states[g,k])
        #         #make sure you didn't end up with any negative numbers due to machine precision
        #         for(j in 1:J.sight[g]){
        #           if(bigLam.marked.proposed[g,j,k]<0){
        #             bigLam.marked.proposed[g,j,k] <- 0
        #           }
        #           if(bigLam.unmarked.proposed[g,j,k]<0){
        #             bigLam.unmarked.proposed[g,j,k] <- 0
        #           }
        #         }
        #       }
        #     }
        #   }
        #   model$calculate(pd.nodes) #update pd nodes
        #   model$calculate(lam.nodes) #update lam nodes
        #   #add these in after calculating lam
        #   for(g2 in 1:n.sight.g){
        #     g <- sight.g[g2]
        #     if(model$z[i,g]==1){
        #       for(k in 1:K.sight[g]){
        #         bigLam.marked.proposed[g,1:J.sight[g],k] <- bigLam.marked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]*mark.states[g,k]
        #         bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]*(1-mark.states[g,k])
        #       }
        #     }
        #   }
        #   #put bigLam.marked in model object
        #   model$bigLam.marked <<- bigLam.marked.proposed
        #   model$bigLam.unmarked <<- bigLam.unmarked.proposed
        #   model$calculate(lam.mnoID.nodes) #update after bigLam
        #   model$calculate(lam.um.nodes) #update after bigLam
        #   model$calculate(lam.unk.nodes) #update after bigLam
        #   lp.proposed.y.mark <- model$calculate(y.mark.nodes)
        #   lp.proposed.y.mID <- model$calculate(y.mID.nodes)
        #   lp.proposed.y.mnoID <- model$calculate(y.mnoID.nodes)
        #   lp.proposed.y.um <- model$calculate(y.um.nodes)
        #   lp.proposed.y.unk <- model$calculate(y.unk.nodes)
        #   lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.mID + lp.initial.y.mnoID + lp.initial.y.um + lp.initial.y.unk
        #   lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.mID + lp.proposed.y.mnoID + lp.proposed.y.um + lp.proposed.y.unk
        # }else{ #else unmarked
        #   #pull this out of model object
        #   bigLam.unmarked.initial <- model$bigLam.unmarked
        #   #update proposed s
        #   model$s[i, 1:2] <<- s.cand
        #   lp.proposed.s <- model$calculate(s.nodes) #proposed logprob for s.nodes
        #   #subtract these out before calculating lam
        #   bigLam.unmarked.proposed <- bigLam.unmarked.initial
        #   for(g2 in 1:n.sight.g){ #z.super always 1 here
        #     g <- sight.g[g2]
        #     if(model$z[i,g]==1){
        #       for(k in 1:K.sight[g]){
        #         bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] - model$lam[i,g,1:J.sight[g]]
        #         #make sure you didn't end up with any negative numbers due to machine precision
        #         for(j in 1:J.sight[g]){
        #           if(bigLam.unmarked.proposed[g,j,k]<0){
        #             bigLam.unmarked.proposed[g,j,k] <- 0
        #           }
        #         }
        #       }
        #     }
        #   }
        #   model$calculate(pd.nodes) #update pd nodes
        #   model$calculate(lam.nodes) #update lam nodes
        #   #add these in after calculating lam
        #   for(g2 in 1:n.sight.g){
        #     g <- sight.g[g2]
        #     if(model$z[i,g]==1){
        #       for(k in 1:K.sight[g]){
        #         bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]
        #       }
        #     }
        #   }
        #   #put bigLam in model object
        #   model$bigLam.unmarked <<- bigLam.unmarked.proposed
        #   model$calculate(lam.um.nodes) #update after bigLam
        #   model$calculate(lam.unk.nodes) #update after bigLam
        #   lp.proposed.y.mark <- model$calculate(y.mark.nodes)
        #   lp.proposed.y.um <- model$calculate(y.um.nodes)
        #   lp.proposed.y.unk <- model$calculate(y.unk.nodes)
        #   lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk
        #   lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk
        # }
        # initial log probability: always-active s/density/telemetry terms
        # plus focal and aggregate observation terms only in primary occasions where z[i,g]=1
        lp.initial <- model$getLogProb(s.nodes)
        for(g in 1:n.primary){
          if(model$z[i,g]==1){
            if(focal.node.n[g]>0){
              node.start <- focal.node.start[g]
              node.end <- node.start+focal.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(focal.nodes[node.start:node.end])
            }
            if(aggregate.node.n[g]>0){
              node.start <- aggregate.node.start[g]
              node.end <- node.start+aggregate.node.n[g]-1
              lp.initial <- lp.initial+model$getLogProb(aggregate.nodes[node.start:node.end])
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
            for(k in 1:K.sight[g]){
              if(mark.states[g,k]==1){
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
              }else{
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
        
        #update proposed s and all focal observation terms in primary occasions where z[i,g]=1
        model$s[i, 1:2] <<- s.cand
        lp.proposed <- model$calculate(s.nodes)
        for(g in 1:n.primary){
          if(model$z[i,g]==1&focal.node.n[g]>0){
            node.start <- focal.node.start[g]
            node.end <- node.start+focal.node.n[g]-1
            lp.proposed <- lp.proposed+model$calculate(focal.nodes[node.start:node.end])
          }
        }
        
        #add this individual's new lam contribution back into the aggregate rates
        for(g2 in 1:n.sight.g){
          g <- sight.g[g2]
          if(model$z[i,g]==1){
            for(k in 1:K.sight[g]){
              bigLam.marked.proposed[g,1:J.sight[g],k] <- bigLam.marked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]*mark.states[g,k]
              bigLam.unmarked.proposed[g,1:J.sight[g],k] <- bigLam.unmarked.proposed[g,1:J.sight[g],k] + model$lam[i,g,1:J.sight[g]]*(1-mark.states[g,k])
            }
          }
        }
        model$bigLam.marked <<- bigLam.marked.proposed
        model$bigLam.unmarked <<- bigLam.unmarked.proposed
        
        #update aggregate observation terms only in primary occasions where this individual's contribution changed
        for(g in 1:n.primary){
          if(model$z[i,g]==1&aggregate.node.n[g]>0){
            node.start <- aggregate.node.start[g]
            node.end <- node.start+aggregate.node.n[g]-1
            lp.proposed <- lp.proposed+model$calculate(aggregate.nodes[node.start:node.end])
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
