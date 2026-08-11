dHabYear1 <- nimbleFunction(
  run = function(x = double(1),pi.cell = double(1),cells = double(2),res = double(0),
                 dSS = double(2),xlim = double(1),ylim = double(1),z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==1){
      # uniform over state space cell
      logProb.unif <- - 2*log(res)
      # cell likelihood
      cell <- cells[trunc(x[1]/res)+1, trunc(x[2]/res)+1]
      logProb.cell <- log(pi.cell[cell])
      logProb <- logProb.unif + logProb.cell
    }else{
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rHabYear1 <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(1),cells = double(2),res = double(0),
                 dSS = double(2),xlim = double(1),ylim = double(1),z.super = double(0)){
    returnType(double(1))
    if(z.super==1){
      #simulate cell from pi.cell then uniform within cell
      n.cells <- nimDim(pi.cell)[1]
      s.cell <- rcat(1, pi.cell[1:n.cells])
      x.sim <- runif(1, dSS[s.cell,1] - res/2, dSS[s.cell,1] + res/2)
      y.sim <- runif(1, dSS[s.cell,2] - res/2, dSS[s.cell,2] + res/2)
      return(c(x.sim, y.sim))
    }else{
      return(c(0,0))
    }
  }
)

dHabMove <- nimbleFunction(
  run = function(x = double(1),s.prev = double(1),use.dist = double(1),dSS = double(2),
                 cells = double(2),res = double(0),sigma.move = double(0),z.super = double(0),
                 log = integer(0)){
    returnType(double(0))
    if(z.super==1){
      #find cell with lookup table
      cell <- cells[trunc(x[1]/res)+1,trunc(x[2]/res)+1]
      #get cell bounds
      x.min <- dSS[cell,1] - res/2
      x.max <- dSS[cell,1] + res/2
      y.min <- dSS[cell,2] - res/2
      y.max <- dSS[cell,2] + res/2
      # cell selection logProb
      logProb.cell <- log(use.dist[cell])
      # continuous within-cell location logProb
      logProb.x <- dnorm(x[1],s.prev[1],sigma.move, log=TRUE) -
        log(pnorm(x.max,s.prev[1],sigma.move) - pnorm(x.min,s.prev[1],sigma.move))
      logProb.y <- dnorm(x[2],s.prev[2],sigma.move,log=TRUE) -
        log(pnorm(y.max,s.prev[2],sigma.move) - pnorm(y.min,s.prev[2],sigma.move))
      logProb <- logProb.cell + logProb.x + logProb.y
    }else{
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rHabMove <- nimbleFunction(
  run = function(n = integer(0),s.prev = double(1),use.dist = double(1),dSS = double(2),
                 cells = double(2),res = double(0),sigma.move = double(0),z.super = double(0)) {
    returnType(double(1))
    if(z.super==1){
      n.cells <- nimDim(dSS)[1]
      s.cell <- rcat(1,use.dist[1:n.cells])
      x.min <- dSS[s.cell,1] - res/2
      x.max <- dSS[s.cell,1] + res/2
      y.min <- dSS[s.cell,2] - res/2
      y.max <- dSS[s.cell,2] + res/2
      F.a <- pnorm(x.min,s.prev[1],sigma.move)
      F.b <- pnorm(x.max,s.prev[1],sigma.move)
      x.sim <- qnorm(runif(1,F.a,F.b),s.prev[1],sigma.move)
      F.a <- pnorm(y.min,s.prev[2],sigma.move)
      F.b <- pnorm(y.max,s.prev[2],sigma.move)
      y.sim <- qnorm(runif(1,F.a,F.b),s.prev[2],sigma.move)
      return(c(x.sim,y.sim))
    }else{
      return(c(0,0))
    }
  }
)


getAvail <- nimbleFunction(
  run = function(s = double(1),sigma=double(0),res=double(0),x.vals=double(1),y.vals=double(1),
                 n.cells.x=integer(0),n.cells.y=integer(0),z.super=double(0)) {
    returnType(double(1))
    if(z.super==1){
      avail.dist.x <- rep(0,n.cells.x)
      avail.dist.y <- rep(0,n.cells.y)
      delta <- 1e-8 #this sets the degree of trimming used to get individual availability distributions
      x.limits <- qnorm(c(delta,1-delta),mean=s[1],sd=sigma)
      y.limits <- qnorm(c(delta,1-delta),mean=s[2],sd=sigma)
      #convert to grid edges instead of centroids
      x.vals.edges <- c(x.vals - res/2, x.vals[n.cells.x]+0.5*res)
      y.vals.edges <- c(y.vals - res/2, y.vals[n.cells.y]+0.5*res)
      #trim in x and y direction
      if(x.vals.edges[1]<x.limits[1]){
        x.start <- floor((x.limits[1] - x.vals.edges[1]) / res) + 1
      }else{
        x.start <- 1
      }
      if(x.vals.edges[n.cells.x]>x.limits[2]){
        x.stop <- ceiling((x.limits[2] - x.vals.edges[1]) / res)
      }else{
        x.stop <- n.cells.x
      }
      if(y.vals.edges[1]<y.limits[1]){
        y.start <- floor((y.limits[1] - y.vals.edges[1]) / res) + 1
      }else{
        y.start <- 1
      }
      if(y.vals.edges[n.cells.y]>y.limits[2]){
        y.stop <- ceiling((y.limits[2] - y.vals.edges[1]) / res)
      }else{
        y.stop <- n.cells.y
      }
      pnorm.x <- rep(0,n.cells.x+1)
      pnorm.y <- rep(0,n.cells.y+1)
      #get pnorms
      for(l in x.start:(x.stop+1)){
        pnorm.x[l] <- pnorm(x.vals.edges[l],mean=s[1],sd=sigma)
      }
      for(l in y.start:(y.stop+1)){
        pnorm.y[l] <- pnorm(y.vals.edges[l],mean=s[2],sd=sigma)
      }
      for(l in (x.start):(x.stop)){
        avail.dist.x[l] <- pnorm.x[l+1] - pnorm.x[l]
      }
      for(l in (y.start):(y.stop)){
        avail.dist.y[l] <- pnorm.y[l+1] - pnorm.y[l]
      }
      avail.dist.tmp <- matrix(0,n.cells.x,n.cells.y)
      sum.dist <- 0
      for(i in x.start:x.stop){
        for(j in y.start:y.stop){
          avail.dist.tmp[i,j] <- avail.dist.x[i]*avail.dist.y[j]
          sum.dist <- sum.dist + avail.dist.tmp[i,j]
        }
      }
      avail.dist <- c(avail.dist.tmp)
      #if any probability mass is outside state space, normalize
      if(sum.dist<1){
        avail.dist <- avail.dist/sum.dist
      }
    }else{
      n.cells <- n.cells.x*n.cells.y
      avail.dist <- rep(0,n.cells)
    }
    return(avail.dist)
  }
)

getUse <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1),z.super = double(0)){
    returnType(double(1))
    if(z.super==1){
      use.dist <- rsf*avail.dist
      use.dist <- use.dist/sum(use.dist)
    }else{
      n.cells <- nimDim(rsf)[1]
      use.dist <- rep(0,n.cells)
    }
    return(use.dist)
  }
)

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
    n.primary <- length(phi)
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
#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    n.marked.all <- control$n.marked.all
    n.cap.all <- control$n.cap.all
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    K.sight <- control$K.sight
    xlim <- control$xlim
    ylim <- control$ylim
    x.vals <- control$x.vals
    y.vals <- control$y.vals
    cells <- control$cells
    dSS <- control$dSS
    n.cells <- control$n.cells
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    res <- control$res
    mark.states <- control$mark.states
    mark.states2D <- control$mark.states2D
    tel.z.states <- control$tel.z.states
    y2D <- control$y2D
    mark.years <- control$mark.years
    sight.years <- control$sight.years
    n.mark.years <- control$n.mark.years
    n.sight.years <- control$n.sight.years
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
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    s.nodes <- control$s.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    entry.counts.curr <- rep(0,n.primary+1)
    for(g in 1:n.primary){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1)
    }
    entry.counts.curr[n.primary + 1] <- sum(model$z.super==0)
    
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
        i.idx.mark <- seq(i,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
        i.idx.sight <- seq(i,M*n.sight.years,M) #used to reference correct sighting process nodes (y.um, y.unk and lam nodes
        
        #pull this out of model object
        bigLam.unmarked.initial <- model$bigLam.unmarked
        #subtract out this individual's lambdas
        bigLam.unmarked.removed <- bigLam.unmarked.initial
        for(g in 1:n.sight.years){
          gg <- sight.years[g]
          if(z.curr[gg]==1&mark.states2D[i,gg]==0){
            for(k in 1:K.sight[gg]){
              if(mark.states[i,gg,k]==0){
                bigLam.unmarked.removed[gg,1:J.sight[gg],k] <- bigLam.unmarked.removed[gg,1:J.sight[gg],k] - model$lam[i,gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  if(bigLam.unmarked.removed[gg,j,k]<0){
                    bigLam.unmarked.removed[gg,j,k] <- 0
                  }
                }
              }
            }
          }
        }
        
        for(g in 1:first.det){ #must be recruited in primary occasion with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes when a z changes
          
          # add in this individual's lambdas for this z.prop
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          for(g2 in 1:n.sight.years){
            gg <- sight.years[g2]
            if(z.prop[gg]==1&mark.states2D[i,gg]==0){
              for(k in 1:K.sight[gg]){
                if(mark.states[i,gg,k]==0){
                  bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] + model$lam[i,gg,1:J.sight[gg]]
                }
              }
            }
          }
          
          #put bigLam.unmarked in model object
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes) #update after bigLam
          model$calculate(lam.unk.nodes) #update after bigLam
          
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.mark.nodes[i.idx.mark]) + 
            model$calculate(y.um.nodes) +
            model$calculate(y.unk.nodes)
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.tel.z.states + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes
          #add in this individual's lambdas for this z.prop
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          for(g2 in 1:n.sight.years){
            gg <- sight.years[g2]
            if(z.prop[gg]==1&mark.states2D[i,gg]==0){
              for(k in 1:K.sight[gg]){
                if(mark.states[i,gg,k]==0){
                  bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] + model$lam[i,gg,1:J.sight[gg]]
                }
              }
            }
          }
          #put bigLam.unmarked in model object
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes) #update after bigLam
          model$calculate(lam.unk.nodes) #update after bigLam
          
          #update these logProbs
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.um.nodes)
          model$calculate(y.unk.nodes)
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]]
            mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]]
            mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]]
            for(j in 1:J.sight[gg]){
              mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
            }
          }
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]]
            model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]]
            model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]]
            for(j in 1:J.sight[gg]){
              model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.um.nodes)
          model$calculate(y.unk.nodes)
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
        i.idx.mark <- seq(i,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
        i.idx.sight <- seq(i,M*n.sight.years,M) #used to reference correct sighting process nodes (y.um, y.unk and lam nodes
        #pull these out of model object
        bigLam.unmarked.initial <- model$bigLam.unmarked
        #subtract out this individual's lambdas
        bigLam.unmarked.removed <- bigLam.unmarked.initial
        for(g in 1:n.sight.years){
          gg <- sight.years[g]
          if(z.curr[gg]==1&mark.states2D[i,gg]==0){
            for(k in 1:K.sight[gg]){
              if(mark.states[i,gg,k]==0){
                bigLam.unmarked.removed[gg,1:J.sight[gg],k] <- bigLam.unmarked.removed[gg,1:J.sight[gg],k] - model$lam[i,gg,1:J.sight[gg]]
                for(j in 1:J.sight[gg]){
                  if(bigLam.unmarked.removed[gg,j,k]<0){
                    bigLam.unmarked.removed[gg,j,k] <- 0
                  }
                }
              }
            }
          }
        }
        for(g in (last.det):n.primary){ #can't die on or before primary occasion of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes when a z changes
          # add in this individual's lambdas for this z.prop
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          for(g2 in 1:n.sight.years){
            gg <- sight.years[g2]
            if(z.prop[gg]==1&mark.states2D[i,gg]==0){
              for(k in 1:K.sight[gg]){
                if(mark.states[i,gg,k]==0){
                  bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] + model$lam[i,gg,1:J.sight[gg]]
                }
              }
            }
          }
          #put bigLam.marked and bigLam.unmarked in model object
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes) #update after bigLam
          model$calculate(lam.unk.nodes) #update after bigLam
          
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.mark.nodes[i.idx.mark]) + 
            model$calculate(y.um.nodes) +
            model$calculate(y.unk.nodes)
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.tel.z.states
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
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes
          #add in this individual's lambdas for this z.prop
          bigLam.unmarked.proposed <- bigLam.unmarked.removed
          for(g2 in 1:n.sight.years){
            gg <- sight.years[g2]
            if(z.prop[gg]==1&mark.states2D[i,gg]==0){
              for(k in 1:K.sight[gg]){
                if(mark.states[i,gg,k]==0){
                  bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] + model$lam[i,gg,1:J.sight[gg]]
                }
              }
            }
          }
          #put bigLam.marked and bigLam.unmarked in model object
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes) #update after bigLam
          model$calculate(lam.unk.nodes) #update after bigLam
          #update these logProbs
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.um.nodes)
          model$calculate(y.unk.nodes)
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]]
            mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]]
            mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]]
            for(j in 1:J.sight[gg]){
              mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]]
            model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]]
            model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]]
            for(j in 1:J.sight[gg]){
              model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.um.nodes)
          model$calculate(y.unk.nodes)
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
        # i.idx <- seq(i,M*n.primary,M) #used to reference correct y and pd nodes
        i.idx.mark <- seq(i,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
        i.idx.sight <- seq(i,M*n.sight.years,M) #used to reference correct sighting process nodes (y.um, y.unk and lam nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #get initial logProbs
        lp.initial.entry <- model$getLogProb(N.nodes[1])
        lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
        lp.initial.y.mark <- model$getLogProb(y.mark.nodes[i.idx.mark])
        lp.initial.y.um <- model$getLogProb(y.um.nodes)
        lp.initial.y.unk <- model$getLogProb(y.unk.nodes)
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[i])
        log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))
        
        #track proposal probs - survival is symmetric, but not recruitment and detection
        log.prop.for <- log.prop.back <- 0
        
        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.primary)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])
        
        #simulate survival
        if(z.start.prop < n.primary){#if you don't recruit in final primary occasion
          for(g in (z.start.prop+1):n.primary){
            z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
            log.prop.for <- log.prop.for + dbinom(z.prop[g],1,model$phi[i,g-1]*z.prop[g-1],log=TRUE)
          }
        }
        z.on.prop <- which(z.prop==1)
        z.stop.prop <- max(z.on.prop)
        model$z[i,] <<- z.prop
        model$z.start[i] <<- z.start.prop
        model$z.stop[i] <<- z.stop.prop
        
        #update N, N.recruit, N.survive only if individual is in superpopulation
        #1) Update N
        model$N <<- model$N - z.curr + z.prop
        #2) Update N.recruit
        if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
          model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
        }
        if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
          model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
        
        model$calculate(ER.nodes) #update ER when N updated
        
        #subtract these out before calculating lam, if z.super=1
        bigLam.unmarked.proposed <- bigLam.unmarked.initial
        for(g in 1:n.sight.years){
          gg <- sight.years[g]
          if(z.curr[gg]==1){
            for(k in 1:K.sight[gg]){
              bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] - model$lam[i,gg,1:J.sight[gg]]
              for(j in 1:J.sight[gg]){
                if(bigLam.unmarked.proposed[gg,j,k]<0){
                  bigLam.unmarked.proposed[gg,j,k] <- 0
                }
              }
            }
          }
        }
        model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
        model$calculate(lam.nodes[i.idx.sight]) #update lam nodes after z changes
        #add these in after calculating lam
        for(g in 1:n.sight.years){
          gg <- sight.years[g]
          if(z.prop[gg]==1){
            for(k in 1:K.sight[gg]){
              bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] + model$lam[i,gg,1:J.sight[gg]]
            }
          }
        }
        model$bigLam.unmarked <<- bigLam.unmarked.proposed #put this back in model
        model$calculate(lam.um.nodes)
        model$calculate(lam.unk.nodes)
        
        #get proposed logProbs
        lp.proposed.entry <- model$calculate(N.nodes[1])
        lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
        lp.proposed.y.mark <- model$calculate(y.mark.nodes[i.idx.mark])
        lp.proposed.y.um <- model$calculate(y.um.nodes)
        lp.proposed.y.unk <- model$calculate(y.unk.nodes)
        lp.proposed.surv <- model$calculate(z.nodes[i])
        lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[i])
        
        # Full multinomial coefficient prior for proposed configuration
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
        entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
        log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
        
        #get backwards proposal probs
        recruit.probs.back <- c(model$lambda.y1,model$ER)
        recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
        log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
        if(z.start.curr < n.primary){#if you don't recruit in final primary occasion
          for(g in (z.start.curr+1):n.primary){
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[i,g-1]*z.curr[g-1],log=TRUE)
          }
        }
        lp.initial.total <- lp.initial.entry + lp.initial.y.mark + lp.initial.y.um +
          lp.initial.y.unk + lp.initial.surv + lp.initial.tel.z.states + log.prior.curr
        lp.proposed.total <- lp.proposed.entry + lp.proposed.y.mark + lp.proposed.y.um +
          lp.proposed.y.unk + lp.proposed.surv + lp.proposed.tel.z.states + log.prior.prop
        
        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        # log_MH_ratio <- (lp.proposed) - (lp.initial)
        accept <- decide(log_MH_ratio)
        
        if(accept){
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]]
            mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]]
            mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]]
            for(j in 1:J.sight[gg]){
              mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
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
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]]
            model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]]
            model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]]
            for(j in 1:J.sight[gg]){
              model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
            }
          }
          #set these logProbs back
          model$calculate(N.recruit.nodes)
          model$calculate(N.nodes[1])
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.um.nodes)
          model$calculate(y.unk.nodes)
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    #3) update z.super: Metropolis-Hastings. only involves unmarked individuals
    #entry counts current coming out of undetected ind update
    bigLam.unmarked.initial <- model$bigLam.unmarked #pull this out.
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z.super==1)
        non.init <- length(z.on)
        pick <- rcat(1,rep(1/non.init,non.init))
        pick <- z.on[pick]
        #prereject turning off any marked individuals or if there is a single unmarked individual
        if(model$N.super[1]==(n.cap.all+1)|pick<=n.cap.all){
          reject <- TRUE
        }
        if(!reject){
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          s.curr <- model$s[pick,1:n.primary,1:2]
          use.dist.curr <- model$use.dist[pick,1:(n.primary-1),1:n.cells]
          
          #p select off guy
          log.p.select.for <- log(1/non.init)
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx.mark <- seq(pick,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
          pick.idx.sight <- seq(pick,M*n.sight.years,M)
          pick.idx.s <- seq(pick,M*n.primary,M)
          
          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes[pick.idx.mark])
          lp.initial.y.um <- model$getLogProb(y.um.nodes)
          lp.initial.y.unk <- model$getLogProb(y.unk.nodes)
          lp.initial.surv <- model$getLogProb(z.nodes[pick])
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])
          lp.initial.s <- model$getLogProb(s.nodes[pick.idx.s])
          
          # propose new N.super/z.super/z.start/z.stop
          model$N.super <<-  model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.primary)
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          for(g in 1:n.primary){
            model$s[pick,g,1:2] <<- c(0,0)
          }
          for(g in 1:(n.primary-1)){
            model$avail.dist[pick,g,] <<- rep(0,n.cells)
            model$use.dist[pick,g,] <<- rep(0,n.cells)
          }
          #subtract these out before calculating lam
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            if(z.curr[gg]==1){
              for(k in 1:K.sight[gg]){
                bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] - model$lam[pick,gg,1:J.sight[gg]]
                #make sure you didn't end up with any negative numbers due to machine precision
                for(j in 1:J.sight[gg]){
                  if(bigLam.unmarked.proposed[gg,j,k]<0){
                    bigLam.unmarked.proposed[gg,j,k] <- 0
                  }
                }
              }
            }
          }
          model$calculate(pd.nodes[pick.idx.mark]) #turn pd off
          model$calculate(lam.nodes[pick.idx.sight]) #turn lam off
          model$bigLam.unmarked <<- bigLam.unmarked.proposed #put this back in model
          model$calculate(lam.um.nodes)
          model$calculate(lam.unk.nodes)
          
          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          if(z.start.curr < n.primary){
            for(g in (z.start.curr+1):n.primary){
              log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
            }
          }
          log.prop.back.s <- dHabYear1(s.curr[1,1:2],pi.cell=model$pi.cell[1:n.cells],
                                       cells=cells[1:n.cells.x,1:n.cells.y],
                                       dSS=dSS[1:n.cells,1:2],res=res,
                                       xlim=xlim,ylim=ylim,z.super=1,log=TRUE)
          for(g in 2:n.primary){
            log.prop.back.s <- log.prop.back.s + dHabMove(s.curr[g,1:2],s.prev=s.curr[g-1,1:2],
                                                          use.dist=use.dist.curr[g-1,1:n.cells],
                                                          dSS=dSS[1:n.cells,1:2],
                                                          cells=cells[1:n.cells.x,1:n.cells.y],
                                                          res=res,sigma.move=model$sigma.move[1],
                                                          z.super=1,log=TRUE)
          }
          
          #get proposed logProbs for N, N.recruit, and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- model$calculate(y.mark.nodes[pick.idx.mark]) #will always be 0
          lp.proposed.y.um <- model$calculate(y.um.nodes)
          lp.proposed.y.unk <- model$calculate(y.unk.nodes)
          lp.proposed.surv <- model$calculate(z.nodes[pick]) #will always be 0
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick]) #will always be 0
          lp.proposed.s <- model$calculate(s.nodes[pick.idx.s]) #will always be 0
          
          lp.initial.total <- lp.initial.N + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk +
            lp.initial.N.recruit + lp.initial.surv + lp.initial.tel.z.states + lp.initial.s
          lp.proposed.total <- lp.proposed.N + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk +
            lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.tel.z.states + lp.proposed.s
          
          #backwards prior and select probs
          #move from class z.start.curr in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.primary + 1] <- entry.counts.prop[n.primary + 1] + 1
          
          #p select on guy
          noff.back <- sum(model$z.super == 0)
          log.p.select.back <- log(1/noff.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.for <- 0
          log.prop.for.s <- 0
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back + log.prop.back.s) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for + log.prop.for.s)
          
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            mvSaved["s",1][pick,1:n.primary,1:2] <<- model[["s"]][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              mvSaved["avail.dist",1][pick,g,] <<- model[["avail.dist"]][pick,g,]
              mvSaved["use.dist",1][pick,g,] <<- model[["use.dist"]][pick,g,]
            }
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]]
              mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]]
              mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
              }
            }
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            model[["s"]][pick,1:n.primary,1:2] <<- mvSaved["s",1][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              model[["avail.dist"]][pick,g,] <<- mvSaved["avail.dist",1][pick,g,]
              model[["use.dist"]][pick,g,] <<- mvSaved["use.dist",1][pick,g,]
            }
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]]
              model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]]
              model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]]
              for(j in 1:J.sight[gg]){
                model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
              }
            }
            #set these logProbs back
            model$calculate(y.mark.nodes[pick.idx.mark])
            model$calculate(y.um.nodes)
            model$calculate(y.unk.nodes)
            model$calculate(z.nodes[pick])
            model$calculate(tel.z.states.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            model$calculate(s.nodes[pick.idx.s])
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          z.off <- which(model$z.super==0)
          noff.init <- length(z.off)
          pick <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx.mark <- seq(pick,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
          pick.idx.sight <- seq(pick,M*n.sight.years,M)
          pick.idx.s <- seq(pick,M*n.primary,M)
          
          non.init <- sum(model$z.super == 1)
          
          #p select off guy
          log.p.select.for <- log(1/noff.init)
          
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          
          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes[pick.idx.mark]) #will always be 0
          lp.initial.y.um <- model$getLogProb(y.um.nodes)
          lp.initial.y.unk <- model$getLogProb(y.unk.nodes) 
          lp.initial.surv <- model$getLogProb(z.nodes[pick]) #will always be 0
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick]) #will always be 0
          lp.initial.s <- model$getLogProb(s.nodes[pick.idx.s]) #will always be 0
          
          # Propose new z.start for the new on individual
          recruit.probs.for <- c(model$lambda.y1,model$ER)
          recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
          z.start.prop <- rcat(1,recruit.probs.for)  # propose entry cohort
          log.prop.for <- log(recruit.probs.for[z.start.prop])
          model$z.start[pick] <<- z.start.prop
          
          # Simulate survival path
          model$z[pick,] <<- 0 # initialize to 0
          model$z[pick,z.start.prop] <<- 1
          if(z.start.prop < n.primary){
            for(g in (z.start.prop+1):n.primary){
              model$z[pick,g] <<- rbinom(1,1,model$phi[pick,g-1]*model$z[pick,g-1])
              log.prop.for <- log.prop.for + dbinom(model$z[pick,g],1,model$phi[pick,g-1]*model$z[pick,g-1],log=TRUE)
            }
          }
          # Update z.stop
          z.on.prop <- which(model$z[pick,] == 1)
          z.stop.prop <- max(z.on.prop)
          model$z.stop[pick] <<- z.stop.prop
          
          #propose new N/z
          model$N.super <<-  model$N.super + 1
          model$z.super[pick] <<- 1
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N + model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){ #if wasn't in pop in primary occasion 1
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          model$s[pick,1,1:2] <<- rHabYear1(1,pi.cell=model$pi.cell[1:n.cells],
                                            cells=cells[1:n.cells.x,1:n.cells.y],
                                            dSS=dSS[1:n.cells,1:2],res=res,
                                            xlim=xlim,ylim=ylim,z.super=1)
          log.prop.for.s <- dHabYear1(model$s[pick,1,1:2],pi.cell=model$pi.cell[1:n.cells],
                                      cells=cells[1:n.cells.x,1:n.cells.y],
                                      dSS=dSS[1:n.cells,1:2],res=res,
                                      xlim=xlim,ylim=ylim,z.super=1,log=TRUE)
          for(g in 2:n.primary){
            model$avail.dist[pick,g-1,1:n.cells] <<- getAvail(s=model$s[pick,g-1,1:2],
                                                              sigma=model$sigma.move[1],res=res,
                                                              x.vals=x.vals,y.vals=y.vals,
                                                              n.cells.x=n.cells.x,n.cells.y=n.cells.y,
                                                              z.super=1)
            model$use.dist[pick,g-1,1:n.cells] <<- getUse(rsf=model$rsf[1:n.cells],
                                                          avail.dist=model$avail.dist[pick,g-1,1:n.cells],
                                                          z.super=1)
            model$s[pick,g,1:2] <<- rHabMove(1,s.prev=model$s[pick,g-1,1:2],
                                             use.dist=model$use.dist[pick,g-1,1:n.cells],
                                             dSS=dSS[1:n.cells,1:2],
                                             cells=cells[1:n.cells.x,1:n.cells.y],
                                             res=res,sigma.move=model$sigma.move[1],z.super=1)
            log.prop.for.s <- log.prop.for.s + dHabMove(model$s[pick,g,1:2],s.prev=model$s[pick,g-1,1:2],
                                                        use.dist=model$use.dist[pick,g-1,1:n.cells],
                                                        dSS=dSS[1:n.cells,1:2],
                                                        cells=cells[1:n.cells.x,1:n.cells.y],
                                                        res=res,sigma.move=model$sigma.move[1],
                                                        z.super=1,log=TRUE)
          }
          model$calculate(pd.nodes[pick.idx.mark]) #turn pd on
          model$calculate(lam.nodes[pick.idx.sight]) #turn lam on
          #add these in after calculating lam
          bigLam.unmarked.proposed <- bigLam.unmarked.initial
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            if(model$z[pick,gg]==1){
              for(k in 1:K.sight[gg]){
                bigLam.unmarked.proposed[gg,1:J.sight[gg],k] <- bigLam.unmarked.proposed[gg,1:J.sight[gg],k] + model$lam[pick,gg,1:J.sight[gg]]
              }
            }
          }
          model$bigLam.unmarked <<- bigLam.unmarked.proposed #put this back in model
          model$calculate(lam.um.nodes)
          model$calculate(lam.unk.nodes)
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- model$calculate(y.mark.nodes[pick.idx.mark])
          lp.proposed.y.um <- model$calculate(y.um.nodes)
          lp.proposed.y.unk <- model$calculate(y.unk.nodes)
          lp.proposed.surv <- model$calculate(z.nodes[pick])
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])
          lp.proposed.s <- model$calculate(s.nodes[pick.idx.s])
          
          lp.initial.total <- lp.initial.N + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk +
            lp.initial.N.recruit + lp.initial.surv + lp.initial.tel.z.states + lp.initial.s
          lp.proposed.total <- lp.proposed.N + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk +
            lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.tel.z.states + lp.proposed.s
          
          #backwards prior and select probs
          #move from class g in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.prop[n.primary + 1] <- entry.counts.prop[n.primary + 1] - 1
          
          #p select on guy
          non.back <- sum(model$z.super == 1)
          log.p.select.back <- log(1/non.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.back <- 0
          log.prop.back.s <- 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back + log.prop.back.s) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for + log.prop.for.s)
          
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            mvSaved["s",1][pick,1:n.primary,1:2] <<- model[["s"]][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              mvSaved["avail.dist",1][pick,g,] <<- model[["avail.dist"]][pick,g,]
              mvSaved["use.dist",1][pick,g,] <<- model[["use.dist"]][pick,g,]
            }
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]]
              mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]]
              mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]] <<- model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
              }
            }
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            model[["s"]][pick,1:n.primary,1:2] <<- mvSaved["s",1][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              model[["avail.dist"]][pick,g,] <<- mvSaved["avail.dist",1][pick,g,]
              model[["use.dist"]][pick,g,] <<- mvSaved["use.dist",1][pick,g,]
            }
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              model[["bigLam.unmarked"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["bigLam.unmarked",1][gg,1:J.sight[gg],1:K.sight[gg]]
              model[["lam.um"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.um",1][gg,1:J.sight[gg],1:K.sight[gg]]
              model[["lam.unk"]][gg,1:J.sight[gg],1:K.sight[gg]] <<- mvSaved["lam.unk",1][gg,1:J.sight[gg],1:K.sight[gg]]
              for(j in 1:J.sight[gg]){
                model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
              }
            }
            #set these logProbs back
            model$calculate(y.mark.nodes[pick.idx.mark])
            model$calculate(y.um.nodes)
            model$calculate(y.unk.nodes)
            model$calculate(z.nodes[pick])
            model$calculate(tel.z.states.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            model$calculate(s.nodes[pick.idx.s])
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
    if(target=="gamma"){
      is.fixed.gamma <- TRUE
      g <- 1
      recruitNodes <- grep("^N.recruit\\[",model$getNodeNames(stochOnly=TRUE),value=TRUE)
      n.recruit <- length(recruitNodes)
    }else{
      is.fixed.gamma <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
      n.recruit <- 1
    }
  },
  run = function(){
    if(is.fixed.gamma){
      count <- 0
      rate <- 0
      for(j in 1:n.recruit){
        count <- count+model$N.recruit[j]
        rate <- rate+model$N[j]
      }
    }else{
      count <- model$N.recruit[g]
      rate <- model$N[g]
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