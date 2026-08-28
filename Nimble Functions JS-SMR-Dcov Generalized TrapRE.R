#telemetry survival vector distribution
dSurvivalTel <- nimbleFunction(
  run = function(x = double(1), z = double(1), z.super = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      for(g in 1:length(z)){
        if(x[g]!=2){ #2 codes NA (not yet collared or uninformatively censored), 1 is alive, 0 is dead
          #z invalid if conflicts with tel.
          if(x[g] != z[g]){
            logProb <- -Inf
          }
        }
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rSurvivalTel <- nimbleFunction(
  run = function(n = integer(0), z = double(1),z.super = double(0)) {
    returnType(double(1))
    return(rep(0,length(z)))
  }
)

#telemetry location vector distribution
dNormVector <- nimbleFunction(
  run = function(x = double(2), s = double(1), sigma = double(0), n.locs.ind = double(0),
                 max.n.tel.locs = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(n.locs.ind>0){
      for(i in 1:n.locs.ind){
        logProb <- logProb + dnorm(x[i,1], mean = s[1], sd = sigma, log = TRUE)
        logProb <- logProb + dnorm(x[i,2], mean = s[2], sd = sigma, log = TRUE)
      }
    }
    return(logProb) 
  }
)

rNormVector <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sigma = double(0), n.locs.ind = double(0),
                 max.n.tel.locs = double(0)) {
    returnType(double(2))
    out <- matrix(0,nrow=max.n.tel.locs,ncol=2)
    return(out)
  }
)

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

GetbigLam <- nimbleFunction(
  run = function(lam = double(2), z = double(1)){ 
    returnType(double(1))
    M <- nimDim(lam)[1]
    J <- nimDim(lam)[2]
    bigLam <- rep(0,J)
    for(i in 1:M){
      if(z[i]==1){
        bigLam <- bigLam + lam[i,]
      }
    }
    return(bigLam)
  }
)

#this is used to restrict likelihood evaluation to only the primary occasions relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.primary <- length(phi)+1
      #extract first and last survival event primary occasions
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.primary){ #if surv.start beyond last primary occasion, no survival events, logProb=0
        if(surv.stop > n.primary){ #but can't survive past n.primary
          surv.stop <- n.primary 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event primary occasions
          logProb <- logProb + dbinom(x[g], size = 1, p = phi[g-1], log = TRUE)
        }
      }
    }
    return(logProb)
  }
)

#make dummy random vector generator to make nimble happy
rSurvival <- nimbleFunction(
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0),z.super = double(0)) {
    returnType(double(1))
    n.primary <- length(phi)+1
    return(rep(0,n.primary))
  }
)

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
    }
    if(z==1){ #otherwise calculate
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0), z.super = double(0),mark.states = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if((z.super*z*mark.states)==0){#skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
      return(0)
    }else{
      logProb <- sum(dpois(x, lambda=lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0), lam = double(1), z = double(0), z.super = double(0),mark.states = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

dDetectorNBCorrection <- nimbleFunction(
  run=function(x=double(0), y.trap.total=double(1), bigLam.marked=double(1), bigLam.unmarked=double(1),
               K1D=double(1), theta=double(0), J=integer(0), log=integer(0)){
    returnType(double(0))
    if(theta<=0){
      if(log){
        return(-Inf)
      }else{
        return(0)
      } 
    }
    logProb <- 0
    for(j in 1:J){
      mu <- K1D[j]*(bigLam.marked[j]+bigLam.unmarked[j])
      count <- y.trap.total[j]
      logProb <- logProb + lgamma(theta+count)-lgamma(theta) + theta*log(theta)+mu-(theta+count)*log(theta+mu)
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    } 
  }
)

#dummy RNG to make nimble happy
rDetectorNBCorrection <- nimbleFunction(
  run=function(n=integer(0), y.trap.total=double(1), bigLam.marked=double(1), bigLam.unmarked=double(1),
               K1D=double(1), theta=double(0), J=integer(0)){
    returnType(double(0))
    return(0)
  }
)

# Function to calculate detection rate, but skip when z=0
GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J))
    }else{
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)
#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K1D = double(1), z = double(0),z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, size = K1D, prob = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0),pd = double(1), K1D = double(1), z = double(0),z.super = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    n.marked.all <- control$n.marked.all
    n.cap.all <- control$n.cap.all
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    mark.states <- control$mark.states
    tel.z.states <- control$tel.z.states
    y2D <- control$y2D
    mark.g <- control$mark.g
    sight.g <- control$sight.g
    n.mark.g <- control$n.mark.g
    n.sight.g <- control$n.sight.g
    z.super.ups <- control$z.super.ups
    n.primary <- control$n.primary
    z.nodes <- control$z.nodes
    tel.z.states.nodes <- control$tel.z.states.nodes
    y.mark.nodes <- control$y.mark.nodes
    y.um.nodes <- control$y.um.nodes
    y.unk.nodes <- control$y.unk.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    lam.um.nodes <- control$lam.um.nodes
    lam.unk.nodes <- control$lam.unk.nodes
    trap.RE.nodes <- control$trap.RE.nodes
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    entry.counts.curr <- rep(0,n.primary+1)
    for(i in 1:M){
      if(model$z.super[i]==1){
        z.start.i <- model$z.start[i]
        entry.counts.curr[z.start.i] <- entry.counts.curr[z.start.i]+1
      }else{
        entry.counts.curr[n.primary+1] <- entry.counts.curr[n.primary+1]+1
      }
    }
    
    #1) Detected guy updates: z.start, z.stop
    #detected/observed individuals include actually marked individuals and captured-but-not-marked individuals
    #y.mID logProb does not change because known to be in population those primary occasions
    #same for bigLam.marked and y.mnoID.
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:n.cap.all){
      if(y2D[i,1]==0){ #skip if known to be alive in 1st primary occasion
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y2D[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.primary)
        #only observation nodes before first detection can change across z.start candidates
        mark.before <- which(mark.g<first.det)
        sight.before <- which(sight.g<first.det)
        
        #pull this out of model object
        bigLam.unmarked.initial <- model$bigLam.unmarked
        #subtract out this individual's lambdas only where z can change
        bigLam.unmarked.removed <- bigLam.unmarked.initial
        if(length(sight.before)>0){
          for(g2 in 1:length(sight.before)){
            idx.g <- sight.before[g2]
            gg <- sight.g[idx.g]
            if(z.curr[gg]==1&mark.states[i,gg]==0){
              for(j in 1:J.sight[gg]){
                bigLam.old <- bigLam.unmarked.removed[gg,j]
                bigLam.unmarked.removed[gg,j] <- bigLam.old-model$lam[i,gg,j]
                #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                if(bigLam.old>0&bigLam.unmarked.removed[gg,j]<=1e-12*bigLam.old){
                  bigLam.unmarked.removed[gg,j] <- 0
                  for(k in 1:M){
                    if(k!=i&mark.states[k,gg]==0){
                      bigLam.unmarked.removed[gg,j] <- bigLam.unmarked.removed[gg,j]+model$lam[k,gg,j]
                    }
                  }
                }
                if(bigLam.unmarked.removed[gg,j]<0){
                  bigLam.unmarked.removed[gg,j] <- 0
                }
              }
            }
          }
        }
        
        #remove focal individual from entry counts once. The candidate-specific part of the
        #multinomial coefficient is then just log(entry.counts.minus[g]+1)
        entry.counts.minus <- entry.counts.curr
        entry.counts.minus[z.start.curr] <- entry.counts.minus[z.start.curr]-1
        #all z.start > 1 candidates have the same N[1], so only calculate that logProb once
        lp.N1.not1 <- 0
        
        for(g in 1:first.det){ #must be recruited in primary occasion with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det<n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr-z.curr+z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr>1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1]-1
          }
          if(z.start.prop>1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1]+1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #only ER nodes before first detection can change across z.start candidates
          for(g2 in 1:(first.det-1)){
            model$calculate(ER.nodes[g2])
          }
          #only focal observation nodes before first detection can change
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
          }
          
          #add in this individual's lambdas for this z.prop only where z can change
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              gg <- sight.g[idx.g]
              if(z.prop[gg]==1&mark.states[i,gg]==0){
                bigLam.unmarked.proposed[gg,1:J.sight[gg]] <- bigLam.unmarked.proposed[gg,1:J.sight[gg]]+
                  model$lam[i,gg,1:J.sight[gg]]
              }
            }
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              model$calculate(lam.um.nodes[idx.g]) #update after bigLam
              model$calculate(lam.unk.nodes[idx.g]) #update after bigLam
            }
          }
          
          #get these logProbs
          #there are only two possible N[1] values: z.start=1 and z.start>1
          if(g==1){
            lp.N1 <- model$calculate(N.nodes[1])
          }else{
            if(g==2){
              lp.N1.not1 <- model$calculate(N.nodes[1])
            }
            lp.N1 <- lp.N1.not1
          }
          #only recruitment likelihoods before first detection can change
          lp.N.recruit <- 0
          for(g2 in 1:(first.det-1)){
            lp.N.recruit <- lp.N.recruit+model$calculate(N.recruit.nodes[g2])
          }
          lp.y.mark <- 0
          lp.y.um <- 0
          lp.y.unk <- 0
          lp.trap.RE <- 0
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              lp.y.mark <- lp.y.mark+model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            #these are marginalized likelihood nodes shared by all individuals, so every
            #affected sighting primary occasion must be recalculated
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              lp.y.um <- lp.y.um+model$calculate(y.um.nodes[idx.g])
              lp.y.unk <- lp.y.unk+model$calculate(y.unk.nodes[idx.g])
              lp.trap.RE <- lp.trap.RE+model$calculate(trap.RE.nodes[idx.g])
            }
          }
          lp.surv <- model$calculate(z.nodes[i])
          #telemetry survival can fix additional z states, so retain this full likelihood
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          #after removing this individual, all multinomial coefficient terms common across
          #candidates cancel, leaving log(n.g+1) for candidate entry class g
          lp.prior <- log(entry.counts.minus[g]+1)
          lp.start[g] <- lp.N1+lp.N.recruit+lp.y.mark+lp.y.um+lp.y.unk+lp.trap.RE+lp.surv+lp.tel.z.states+lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[z.start.prop:first.det] <- 1 #must be alive until first detection
          if(first.det<n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr-z.curr+z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr>1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1]-1
          }
          if(z.start.prop>1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1]+1
          }
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          for(g2 in 1:(first.det-1)){
            model$calculate(ER.nodes[g2])
          }
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
          }
          #add in this individual's lambdas for accepted z.start
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              gg <- sight.g[idx.g]
              if(z.prop[gg]==1&mark.states[i,gg]==0){
                bigLam.unmarked.proposed[gg,1:J.sight[gg]] <- bigLam.unmarked.proposed[gg,1:J.sight[gg]]+
                  model$lam[i,gg,1:J.sight[gg]]
              }
            }
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              model$calculate(lam.um.nodes[idx.g])
              model$calculate(lam.unk.nodes[idx.g])
            }
          }
          #update these logProbs
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              model$calculate(y.um.nodes[idx.g])
              model$calculate(y.unk.nodes[idx.g])
              model$calculate(trap.RE.nodes[idx.g])
            }
          }
          model$calculate(N.nodes[1])
          for(g2 in 1:(first.det-1)){
            model$calculate(N.recruit.nodes[g2])
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              gg <- mark.g[idx.g]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
              }
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              gg <- sight.g[idx.g]
              mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg]]
              mvSaved["lam.um",1][gg,1:J.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg]]
              mvSaved["lam.unk",1][gg,1:J.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg]]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
              }
            }
          }
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr]-1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop]+1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              gg <- mark.g[idx.g]
              for(j in 1:J.mark[gg]){
                model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
              }
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              gg <- sight.g[idx.g]
              model[["bigLam.unmarked"]][gg,1:J.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]]
              model[["lam.um"]][gg,1:J.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg]]
              model[["lam.unk"]][gg,1:J.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg]]
              for(j in 1:J.sight[gg]){
                model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
              }
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          for(g2 in 1:(first.det-1)){
            model$calculate(N.recruit.nodes[g2])
          }
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              idx.g <- mark.before[g2]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              idx.g <- sight.before[g2]
              model$calculate(y.um.nodes[idx.g])
              model$calculate(y.unk.nodes[idx.g])
              model$calculate(trap.RE.nodes[idx.g])
            }
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    
    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    #detected/observed individuals include actually marked individuals and captured-but-not-marked individuals
    for(i in 1:n.cap.all){
      if(y2D[i,n.primary]==0){ #skip if known to be alive in final primary occasion
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y2D[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.primary)
        #only observation nodes after last detection can change across z.stop candidates
        mark.after <- which(mark.g>last.det)
        sight.after <- which(sight.g>last.det)
        
        #pull this out of model object
        bigLam.unmarked.initial <- model$bigLam.unmarked
        #subtract out this individual's lambdas only where z can change
        bigLam.unmarked.removed <- bigLam.unmarked.initial
        if(length(sight.after)>0){
          for(g2 in 1:length(sight.after)){
            idx.g <- sight.after[g2]
            gg <- sight.g[idx.g]
            if(z.curr[gg]==1&mark.states[i,gg]==0){
              for(j in 1:J.sight[gg]){
                bigLam.old <- bigLam.unmarked.removed[gg,j]
                bigLam.unmarked.removed[gg,j] <- bigLam.old-model$lam[i,gg,j]
                #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                if(bigLam.old>0&bigLam.unmarked.removed[gg,j]<=1e-12*bigLam.old){
                  bigLam.unmarked.removed[gg,j] <- 0
                  for(k in 1:M){
                    if(k!=i&mark.states[k,gg]==0){
                      bigLam.unmarked.removed[gg,j] <- bigLam.unmarked.removed[gg,j]+model$lam[k,gg,j]
                    }
                  }
                }
                if(bigLam.unmarked.removed[gg,j]<0){
                  bigLam.unmarked.removed[gg,j] <- 0
                }
              }
            }
          }
        }
        
        for(g in last.det:n.primary){ #can't die on or before primary occasion of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:last.det] <- z.curr[1:last.det] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr-z.curr+z.prop
          #only ER nodes after last detection can change across z.stop candidates
          if(last.det<n.primary-1){
            for(g2 in (last.det+1):(n.primary-1)){
              model$calculate(ER.nodes[g2])
            }
          }
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
          }
          #add in this individual's lambdas for this z.prop only where z can change
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              gg <- sight.g[idx.g]
              if(z.prop[gg]==1&mark.states[i,gg]==0){
                bigLam.unmarked.proposed[gg,1:J.sight[gg]] <- bigLam.unmarked.proposed[gg,1:J.sight[gg]]+
                  model$lam[i,gg,1:J.sight[gg]]
              }
            }
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              model$calculate(lam.um.nodes[idx.g])
              model$calculate(lam.unk.nodes[idx.g])
            }
          }
          
          #get these logProbs
          lp.N.recruit <- 0
          if(last.det<n.primary-1){
            for(g2 in (last.det+1):(n.primary-1)){
              lp.N.recruit <- lp.N.recruit+model$calculate(N.recruit.nodes[g2])
            }
          }
          lp.y.mark <- 0
          lp.y.um <- 0
          lp.y.unk <- 0
          lp.trap.RE <- 0
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              lp.y.mark <- lp.y.mark+model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            #these are marginalized likelihood nodes shared by all individuals, so every
            #affected sighting primary occasion must be recalculated
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              lp.y.um <- lp.y.um+model$calculate(y.um.nodes[idx.g])
              lp.y.unk <- lp.y.unk+model$calculate(y.unk.nodes[idx.g])
              lp.trap.RE <- lp.trap.RE+model$calculate(trap.RE.nodes[idx.g])
            }
          }
          lp.surv <- model$calculate(z.nodes[i])
          #telemetry survival can fix z states after the last live observation, so this likelihood is essential here
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N.recruit+lp.y.mark+lp.y.um+lp.y.unk+lp.trap.RE+lp.surv+lp.tel.z.states
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.primary)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:last.det] <- z.curr[1:last.det] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr-z.curr+z.prop
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          if(last.det<n.primary-1){
            for(g2 in (last.det+1):(n.primary-1)){
              model$calculate(ER.nodes[g2])
            }
          }
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
          }
          #add in this individual's lambdas for accepted z.stop
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              gg <- sight.g[idx.g]
              if(z.prop[gg]==1&mark.states[i,gg]==0){
                bigLam.unmarked.proposed[gg,1:J.sight[gg]] <- bigLam.unmarked.proposed[gg,1:J.sight[gg]]+
                  model$lam[i,gg,1:J.sight[gg]]
              }
            }
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              model$calculate(lam.um.nodes[idx.g])
              model$calculate(lam.unk.nodes[idx.g])
            }
          }
          #update these logProbs
          if(last.det<n.primary-1){
            for(g2 in (last.det+1):(n.primary-1)){
              model$calculate(N.recruit.nodes[g2])
            }
          }
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              model$calculate(y.um.nodes[idx.g])
              model$calculate(y.unk.nodes[idx.g])
              model$calculate(trap.RE.nodes[idx.g])
            }
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              gg <- mark.g[idx.g]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
              }
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              gg <- sight.g[idx.g]
              mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg]]
              mvSaved["lam.um",1][gg,1:J.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg]]
              mvSaved["lam.unk",1][gg,1:J.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg]]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
              }
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              gg <- mark.g[idx.g]
              for(j in 1:J.mark[gg]){
                model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
              }
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              gg <- sight.g[idx.g]
              model[["bigLam.unmarked"]][gg,1:J.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]]
              model[["lam.um"]][gg,1:J.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg]]
              model[["lam.unk"]][gg,1:J.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg]]
              for(j in 1:J.sight[gg]){
                model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
              }
            }
          }
          #set these logProbs back
          if(last.det<n.primary-1){
            for(g2 in (last.det+1):(n.primary-1)){
              model$calculate(N.recruit.nodes[g2])
            }
          }
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              idx.g <- mark.after[g2]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              idx.g <- sight.after[g2]
              model$calculate(y.um.nodes[idx.g])
              model$calculate(y.unk.nodes[idx.g])
              model$calculate(trap.RE.nodes[idx.g])
            }
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    
    #2) undetected guy update. Only if in the superpopulation. Must be unmarked guys
    # Metropolis-Hastings, Propose z vectors from priors
    #entry counts current after z.start update
    bigLam.unmarked.initial <- model$bigLam.unmarked #pull this out.
    for(i in (n.cap.all+1):M){
      if(model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #track proposal probs
        #survival proposal probabilities cancel exactly with the survival likelihood because
        #the survival history is proposed from the same survival model used in the target
        log.prop.for <- log.prop.back <- 0
        
        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.primary)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for+log(recruit.probs.for[z.start.prop])
        
        #simulate survival
        #once the individual dies, remaining z's are already 0 so no more simulation is needed
        z.stop.prop <- z.start.prop
        if(z.start.prop<n.primary){#if you don't recruit in final primary occasion
          for(g in (z.start.prop+1):n.primary){
            if(z.prop[g-1]==1){
              z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
              if(z.prop[g]==1){
                z.stop.prop <- g
              }
            }
          }
        }
        
        #if the proposed history is the current history, there is nothing to calculate or update
        if(z.start.prop!=z.start.curr|z.stop.prop!=z.stop.curr){
          #get initial logProbs
          lp.initial.entry <- model$getLogProb(N.nodes[1])+model$getLogProb(N.recruit.nodes)
          #only observation likelihoods in sampled primary occasions where z changed
          lp.initial.y.mark <- 0
          for(g2 in 1:n.mark.g){
            gg <- mark.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              lp.initial.y.mark <- lp.initial.y.mark+model$getLogProb(y.mark.nodes[i+(g2-1)*M])
            }
          }
          lp.initial.y.um <- 0
          lp.initial.y.unk <- 0
          lp.initial.trap.RE <- 0
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              lp.initial.y.um <- lp.initial.y.um+model$getLogProb(y.um.nodes[g2])
              lp.initial.y.unk <- lp.initial.y.unk+model$getLogProb(y.unk.nodes[g2])
              lp.initial.trap.RE <- lp.initial.trap.RE+model$getLogProb(trap.RE.nodes[g2])
            }
          }
          #demographic survival likelihood cancels exactly with backwards survival proposal probability
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[i])
          
          #subtract out this individual's current lambdas
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.prop[gg]!=z.curr[gg]&z.curr[gg]==1){
              for(j in 1:J.sight[gg]){
                bigLam.old <- bigLam.unmarked.proposed[gg,j]
                bigLam.unmarked.proposed[gg,j] <- bigLam.old-model$lam[i,gg,j]
                #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                if(bigLam.old>0&bigLam.unmarked.proposed[gg,j]<=1e-12*bigLam.old){
                  bigLam.unmarked.proposed[gg,j] <- 0
                  for(k in 1:M){
                    if(k!=i&mark.states[k,gg]==0){
                      bigLam.unmarked.proposed[gg,j] <- bigLam.unmarked.proposed[gg,j]+model$lam[k,gg,j]
                    }
                  }
                }
                if(bigLam.unmarked.proposed[gg,j]<0){
                  bigLam.unmarked.proposed[gg,j] <- 0
                }
              }
            }
          }
          
          model$z[i,] <<- z.prop
          model$z.start[i] <<- z.start.prop
          model$z.stop[i] <<- z.stop.prop
          
          #update N, N.recruit, N.survive only if individual is in superpopulation
          #1) Update N
          model$N <<- model$N-z.curr+z.prop
          #2) Update N.recruit
          if(z.start.curr>1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1]-1
          }
          if(z.start.prop>1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1]+1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated

          #only focal observation nodes in sampled primary occasions where z changed
          for(g2 in 1:n.mark.g){
            gg <- mark.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              model$calculate(pd.nodes[i+(g2-1)*M])
            }
          }
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              model$calculate(lam.nodes[i+(g2-1)*M])
            }
          }
          #add in this individual's proposed lambdas
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            # if(z.prop[gg]==1){
            if(z.prop[gg]!=z.curr[gg]&z.prop[gg]==1){
              bigLam.unmarked.proposed[gg,1:J.sight[gg]] <- bigLam.unmarked.proposed[gg,1:J.sight[gg]]+
                model$lam[i,gg,1:J.sight[gg]]
            }
          }
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          #only aggregate rates in sighting primary occasions where z changed
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              model$calculate(lam.um.nodes[g2])
              model$calculate(lam.unk.nodes[g2])
            }
          }
          
          #get proposed logProbs
          lp.proposed.entry <- model$calculate(N.nodes[1])+model$calculate(N.recruit.nodes)
          #only observation likelihoods in sampled primary occasions where z changed
          lp.proposed.y.mark <- 0
          for(g2 in 1:n.mark.g){
            gg <- mark.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              lp.proposed.y.mark <- lp.proposed.y.mark+model$calculate(y.mark.nodes[i+(g2-1)*M])
            }
          }
          lp.proposed.y.um <- 0
          lp.proposed.y.unk <- 0
          lp.proposed.trap.RE <- 0
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.prop[gg]!=z.curr[gg]){
              lp.proposed.y.um <- lp.proposed.y.um+model$calculate(y.um.nodes[g2])
              lp.proposed.y.unk <- lp.proposed.y.unk+model$calculate(y.unk.nodes[g2])
              lp.proposed.trap.RE <- lp.proposed.trap.RE+model$calculate(trap.RE.nodes[g2])
            }
          }
          #telemetry survival does not cancel against the demographic survival proposal
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          
          #local multinomial coefficient ratio
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr]-1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop]+1
          if(z.start.prop!=z.start.curr){
            log.prior.ratio <- log(entry.counts.curr[z.start.prop]+1)-log(entry.counts.curr[z.start.curr])
          }else{
            log.prior.ratio <- 0
          }
          
          #get backwards proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log.prop.back+log(recruit.probs.back[z.start.curr])
          #survival proposal probabilities are not calculated because they cancel exactly
          #with the survival likelihood ratio in the MH ratio
          lp.initial.total <- lp.initial.entry+lp.initial.y.mark+lp.initial.y.um+
            lp.initial.y.unk+lp.initial.trap.RE+lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.entry+lp.proposed.y.mark+lp.proposed.y.um+
            lp.proposed.y.unk+lp.proposed.trap.RE+lp.proposed.tel.z.states
          
          #MH step
          log_MH_ratio <- (lp.proposed.total+log.prior.ratio+log.prop.back)-(lp.initial.total+log.prop.for)
          accept <- decide(log_MH_ratio)
          
          if(accept){
            #update survival logProb once for accepted history; it was not needed when evaluating MH ratio
            model$calculate(z.nodes[i])
            mvSaved["z.start",1][i] <<- model[["z.start"]][i]
            mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
            mvSaved["z",1][i,] <<- model[["z"]][i,]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["ER",1] <<- model[["ER"]]
            #only changed observation nodes need to be saved
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              if(z.prop[gg]!=z.curr[gg]){
                for(j in 1:J.mark[gg]){
                  mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
                }
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(z.prop[gg]!=z.curr[gg]){
                mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg]]
                mvSaved["lam.um",1][gg,1:J.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg]]
                mvSaved["lam.unk",1][gg,1:J.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
                }
              }
            }
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][i] <<- mvSaved["z.start",1][i]
            model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
            model[["z"]][i,] <<- mvSaved["z",1][i,]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["ER"]] <<- mvSaved["ER",1]
            #set these logProbs back
            model$calculate(N.recruit.nodes)
            model$calculate(N.nodes[1])
            #restore only changed observation nodes and their logProbs
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              if(z.prop[gg]!=z.curr[gg]){
                for(j in 1:J.mark[gg]){
                  model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
                }
                model$calculate(y.mark.nodes[i+(g2-1)*M])
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(z.prop[gg]!=z.curr[gg]){
                model[["bigLam.unmarked"]][gg,1:J.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]]
                model[["lam.um"]][gg,1:J.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg]]
                model[["lam.unk"]][gg,1:J.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
                }
                model$calculate(y.um.nodes[g2])
                model$calculate(y.unk.nodes[g2])
                model$calculate(trap.RE.nodes[g2])
              }
            }
            #model$calculate(z.nodes[i]) #not needed because survival logProb was never recalculated for the proposal
            model$calculate(tel.z.states.nodes[i])
          }
        }
      }
    }
    
    #3) update z.super: Metropolis-Hastings. only involves unmarked individuals
    #entry counts current coming out of undetected ind update
    bigLam.unmarked.initial <- model$bigLam.unmarked #pull this out.
    #make lists of currently on/off unmarked guys once, then update after accepted proposals
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in (n.cap.all+1):M){
      if(model$z.super[i]==1){
        non.curr <- non.curr+1
        z.on[non.curr] <- i
      }else{
        noff.curr <- noff.curr+1
        z.off[noff.curr] <- i
      }
    }
    
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        #keep at least one unmarked individual in the superpopulation
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          
          #p select on guy
          log.p.select.for <- log(1/non.init)
          #get initial logProbs
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          #only need to consider when z.curr=1 since those are the only values that can change with this proposal
          lp.initial.y.mark <- 0
          for(g2 in 1:n.mark.g){
            gg <- mark.g[g2]
            if(z.curr[gg]==1){
              lp.initial.y.mark <- lp.initial.y.mark+model$getLogProb(y.mark.nodes[pick+(g2-1)*M])
            }
          }
          lp.initial.y.um <- 0
          lp.initial.y.unk <- 0
          lp.initial.trap.RE <- 0
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.curr[gg]==1){
              lp.initial.y.um <- lp.initial.y.um+model$getLogProb(y.um.nodes[g2])
              lp.initial.y.unk <- lp.initial.y.unk+model$getLogProb(y.unk.nodes[g2])
              lp.initial.trap.RE <- lp.initial.trap.RE+model$getLogProb(trap.RE.nodes[g2])
            }
          }
          #demographic survival likelihood cancels exactly with reverse survival proposal probability
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])
          
          #subtract out this individual's current unmarked lambdas
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.curr[gg]==1){
              for(j in 1:J.sight[gg]){
                bigLam.old <- bigLam.unmarked.proposed[gg,j]
                bigLam.unmarked.proposed[gg,j] <- bigLam.old-model$lam[pick,gg,j]
                #if subtraction nearly cancels the total, recompute the residual to avoid numerical loss
                if(bigLam.old>0&bigLam.unmarked.proposed[gg,j]<=1e-12*bigLam.old){
                  bigLam.unmarked.proposed[gg,j] <- 0
                  for(k in 1:M){
                    if(k!=pick&mark.states[k,gg]==0){
                      bigLam.unmarked.proposed[gg,j] <- bigLam.unmarked.proposed[gg,j]+model$lam[k,gg,j]
                    }
                  }
                }
                if(bigLam.unmarked.proposed[gg,j]<0){
                  bigLam.unmarked.proposed[gg,j] <- 0
                }
              }
            }
          }
          
          # propose new N.super/z.super/z.start/z.stop
          model$N.super <<- model$N.super-1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.primary)
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N-z.curr
          #2) Update N.recruit
          if(z.start.curr>1){ #if wasn't in pop in primary occasion 1
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1]-1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          #only need to consider when z.curr=1 since those are the only values that can change with this proposal
          #When z.super is proposed off, focal pd and lam are known to be zero.
          #Delay synchronizing those deterministic nodes until the proposal is accepted.
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.curr[gg]==1){
              model$calculate(lam.um.nodes[g2])
              model$calculate(lam.unk.nodes[g2])
            }
          }
          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          #survival proposal probability cancels exactly with current survival likelihood
          
          #get proposed logProbs for N, N.recruit, and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          #only need to consider when z.curr=1 since those are the only values that can change with this proposal
          lp.proposed.y.mark <- 0 #all marking likelihood terms are zero when z.super=0
          lp.proposed.y.um <- 0
          lp.proposed.y.unk <- 0
          lp.proposed.trap.RE <- 0
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(z.curr[gg]==1){
              lp.proposed.y.um <- lp.proposed.y.um+model$calculate(y.um.nodes[g2])
              lp.proposed.y.unk <- lp.proposed.y.unk+model$calculate(y.unk.nodes[g2])
              lp.proposed.trap.RE <- lp.proposed.trap.RE+model$calculate(trap.RE.nodes[g2])
            }
          }
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])
          
          #survival target/proposal terms cancel exactly, so they are omitted from the MH totals
          lp.initial.total <- lp.initial.N+lp.initial.y.mark+lp.initial.y.um+lp.initial.y.unk+
            lp.initial.trap.RE+lp.initial.N.recruit+lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N+lp.proposed.y.mark+lp.proposed.y.um+lp.proposed.y.unk+
            lp.proposed.trap.RE+lp.proposed.N.recruit+lp.proposed.tel.z.states
          
          #backwards prior and select probs
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr]-1
          entry.counts.prop[n.primary+1] <- entry.counts.prop[n.primary+1]+1
          noff.back <- noff.curr+1
          log.p.select.back <- log(1/noff.back)
          #local multinomial coefficient ratio for moving from entry class to off class
          log.z.prior.ratio <- log(entry.counts.curr[n.primary+1]+1)-log(entry.counts.curr[z.start.curr])
          log.prop.for <- 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total+log.z.prior.ratio+log.p.select.back+log.prop.back)-
            (lp.initial.total+log.p.select.for+log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #survival logProb was omitted from the MH calculation because it cancels with the proposal;
            #calculate it once now to synchronize the accepted model state
            model$calculate(z.nodes[pick])
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            #only need to consider when z.curr=1 since those are the only values that can change with this proposal
            #synchronize focal marking nodes only where they changed from on to off
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              if(z.curr[gg]==1){
                model$calculate(pd.nodes[pick+(g2-1)*M])
                model$calculate(y.mark.nodes[pick+(g2-1)*M])
                for(j in 1:J.mark[gg]){
                  mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
                }
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(z.curr[gg]==1){
                model$calculate(lam.nodes[pick+(g2-1)*M])
                mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg]]
                mvSaved["lam.um",1][gg,1:J.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg]]
                mvSaved["lam.unk",1][gg,1:J.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
                }
              }
            }
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
            entry.counts.curr <- entry.counts.prop
            #move guy from on list to off list
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            #only need to consider when z.curr=1 since those are the only values that can change with this proposal
            #Focal pd/lam nodes were not recalculated before this rejected remove proposal, so they remain valid.
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(z.curr[gg]==1){
                model[["bigLam.unmarked"]][gg,1:J.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]]
                model[["lam.um"]][gg,1:J.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg]]
                model[["lam.unk"]][gg,1:J.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg]]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(z.curr[gg]==1){
                model$calculate(y.um.nodes[g2])
                model$calculate(y.unk.nodes[g2])
                model$calculate(trap.RE.nodes[g2])
              }
            }
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            #model$calculate(z.nodes[pick]) #not needed because survival logProb was never recalculated for the proposal
            model$calculate(tel.z.states.nodes[pick])
          }
        }
      }else{#add
        noff.init <- noff.curr
        if(noff.init>0){
          pick.pos <- rcat(1,rep(1/noff.init,noff.init))
          pick <- z.off[pick.pos]
          
          #p select off guy
          log.p.select.for <- log(1/noff.init)
          #get initial logProbs
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          #only need to consider when z.prop=1 since those are the only values that can change with this proposal
          lp.initial.y.mark <- 0
          lp.initial.y.um <- 0
          lp.initial.y.unk <- 0
          lp.initial.trap.RE <- 0
          #demographic survival likelihood cancels exactly with forward survival proposal probability
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])
          
          # Propose new z.start for the new on individual
          recruit.probs.for <- c(model$lambda.y1,model$ER)
          recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
          z.start.prop <- rcat(1,recruit.probs.for)  # propose entry cohort
          log.prop.for <- log(recruit.probs.for[z.start.prop])
          model$z.start[pick] <<- z.start.prop
          
          #Simulate survival path
          model$z[pick,] <<- 0 # initialize to 0
          model$z[pick,z.start.prop] <<- 1
          z.stop.prop <- z.start.prop
          if(z.start.prop<n.primary){
            for(g in (z.start.prop+1):n.primary){
              if(model$z[pick,g-1]==1){
                model$z[pick,g] <<- rbinom(1,1,model$phi[pick,g-1]*model$z[pick,g-1])
                if(model$z[pick,g]==1){
                  z.stop.prop <- g
                }
              }
            }
          }
          model$z.stop[pick] <<- z.stop.prop
          
          #get current aggregate sighting likelihood only in proposed alive primary occasions;
          #these are the only aggregate observation terms that can change if the individual is added
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(model$z[pick,gg]==1){
              lp.initial.y.um <- lp.initial.y.um+model$getLogProb(y.um.nodes[g2])
              lp.initial.y.unk <- lp.initial.y.unk+model$getLogProb(y.unk.nodes[g2])
              lp.initial.trap.RE <- lp.initial.trap.RE+model$getLogProb(trap.RE.nodes[g2])
            }
          }
          
          #propose new N/z
          model$N.super <<- model$N.super+1
          model$z.super[pick] <<- 1
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N+model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick]>1){ #if wasn't in pop in primary occasion 1
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1]+1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          #only need to consider when z.prop=1 since those are the only values that can change with this proposal
          for(g2 in 1:n.mark.g){
            gg <- mark.g[g2]
            if(model$z[pick,gg]==1){
              model$calculate(pd.nodes[pick+(g2-1)*M])
            }
          }
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(model$z[pick,gg]==1){
              model$calculate(lam.nodes[pick+(g2-1)*M])
            }
          }
          #add these in after calculating lam
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(model$z[pick,gg]==1){
              bigLam.unmarked.proposed[gg,1:J.sight[gg]] <- bigLam.unmarked.proposed[gg,1:J.sight[gg]]+
                model$lam[pick,gg,1:J.sight[gg]]
            }
          }
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(model$z[pick,gg]==1){
              model$calculate(lam.um.nodes[g2])
              model$calculate(lam.unk.nodes[g2])
            }
          }
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          #only need to consider when z.prop=1 since those are the only values that can change with this proposal
          lp.proposed.y.mark <- 0
          for(g2 in 1:n.mark.g){
            gg <- mark.g[g2]
            if(model$z[pick,gg]==1){
              lp.proposed.y.mark <- lp.proposed.y.mark+model$calculate(y.mark.nodes[pick+(g2-1)*M])
            }
          }
          lp.proposed.y.um <- 0
          lp.proposed.y.unk <- 0
          lp.proposed.trap.RE <- 0
          for(g2 in 1:n.sight.g){
            gg <- sight.g[g2]
            if(model$z[pick,gg]==1){
              lp.proposed.y.um <- lp.proposed.y.um+model$calculate(y.um.nodes[g2])
              lp.proposed.y.unk <- lp.proposed.y.unk+model$calculate(y.unk.nodes[g2])
              lp.proposed.trap.RE <- lp.proposed.trap.RE+model$calculate(trap.RE.nodes[g2])
            }
          }
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])
          
          #survival target/proposal terms cancel exactly, so they are omitted from the MH totals
          lp.initial.total <- lp.initial.N+lp.initial.y.mark+lp.initial.y.um+lp.initial.y.unk+
            lp.initial.trap.RE+lp.initial.N.recruit+lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N+lp.proposed.y.mark+lp.proposed.y.um+lp.proposed.y.unk+
            lp.proposed.trap.RE+lp.proposed.N.recruit+lp.proposed.tel.z.states
          
          #backwards prior and select probs
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop]+1
          entry.counts.prop[n.primary+1] <- entry.counts.prop[n.primary+1]-1
          non.back <- non.curr+1
          log.p.select.back <- log(1/non.back)
          #local multinomial coefficient ratio for moving from off class to entry class
          log.z.prior.ratio <- log(entry.counts.curr[z.start.prop]+1)-log(entry.counts.curr[n.primary+1])
          log.prop.back <- 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total+log.z.prior.ratio+log.p.select.back+log.prop.back)-
            (lp.initial.total+log.p.select.for+log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #survival logProb was omitted from the MH calculation because it cancels with the proposal;
            #calculate it once now to synchronize the accepted model state
            model$calculate(z.nodes[pick])
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            #only need to consider when z.prop=1 since those are the only values that can change with this proposal
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              if(model$z[pick,gg]==1){
                for(j in 1:J.mark[gg]){
                  mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
                }
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(model$z[pick,gg]==1){
                mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg]]
                mvSaved["lam.um",1][gg,1:J.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg]]
                mvSaved["lam.unk",1][gg,1:J.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
                }
              }
            }
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
            entry.counts.curr <- entry.counts.prop
            #move guy from off list to on list
            z.off[pick.pos] <- z.off[noff.curr]
            z.off[noff.curr] <- 0
            noff.curr <- noff.curr-1
            non.curr <- non.curr+1
            z.on[non.curr] <- pick
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            #only need to consider when z.prop=1 since those are the only values that can change with this proposal
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              if(gg>=z.start.prop&gg<=z.stop.prop){
                for(j in 1:J.mark[gg]){
                  model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
                }
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(gg>=z.start.prop&gg<=z.stop.prop){
                model[["bigLam.unmarked"]][gg,1:J.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg]]
                model[["lam.um"]][gg,1:J.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg]]
                model[["lam.unk"]][gg,1:J.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
                }
              }
            }
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              if(gg>=z.start.prop&gg<=z.stop.prop){
                model$calculate(y.mark.nodes[pick+(g2-1)*M])
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              if(gg>=z.start.prop&gg<=z.stop.prop){
                model$calculate(y.um.nodes[g2])
                model$calculate(y.unk.nodes[g2])
                model$calculate(trap.RE.nodes[g2])
              }
            }
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            #model$calculate(z.nodes[pick]) #not needed because survival logProb was never recalculated for the proposal
            model$calculate(tel.z.states.nodes[pick])
          }
        }
      }
    }
    
    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

truncGammaPoisSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    calcNodes <- model$getDependencies(target)
    upper <- model$getBound(target,"upper")
    tau <- control$tau
    if(target=="lambda.y1"){
      is.lambda <- TRUE
      is.fixed.gamma <- FALSE
      g <- 1
      n.recruit <- 1
    }else if(target=="gamma"){
      is.lambda <- FALSE
      is.fixed.gamma <- TRUE
      g <- 1
      recruitNodes <- grep("^N.recruit\\[",model$getNodeNames(stochOnly=TRUE),value=TRUE)
      n.recruit <- length(recruitNodes)
    }else{
      is.lambda <- FALSE
      is.fixed.gamma <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
      n.recruit <- 1
    }
  },
  run = function(){
    if(is.lambda){
      count <- model$N[1]
      rate <- 1
    }else if(is.fixed.gamma){
      count <- 0
      rate <- 0
      for(j in 1:n.recruit){
        count <- count+model$N.recruit[j]
        rate <- rate+model$N[j]*tau[j]
      }
    }else{
      count <- model$N.recruit[g]
      rate <- model$N[g]*tau[g]
    }
    if(rate>0){
      shape <- count+1
      p.upper <- pgamma(upper,shape=shape,rate=rate)
      model[[target]] <<- qgamma(runif(1,0,p.upper),shape=shape,rate=rate)
    }else{
      model[[target]] <<- runif(1,0,upper)
    }
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)
