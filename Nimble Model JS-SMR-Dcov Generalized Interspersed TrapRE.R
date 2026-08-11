NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  # D.intercept <- exp(D.beta0)*cellArea
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) #separate this component so s's do not depend on D.intercept
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  ##Abundance##
  lambda.y1 <- D.intercept*pi.denom #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.primary){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #abundance by primary occasion
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.primary-1)]) #size of superpopulation
  
   #Recruitment
  gamma ~ dunif(0,2) #fixed recruitment parameter
  for(g in 1:(n.primary-1)){
    ER[g] <- N[g]*gamma #expected recruits, if gamma fixed
    # gamma[g] ~ dunif(0,2) #recruitment priors by primary occasion
    # ER[g] <- N[g]*gamma[g] #expected recruits, if gamma 
    N.recruit[g] ~ dpois(ER[g]) #realized recruits
  }
  
  #Individual covariates
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]]) #categorical likelihood for this cell, equivalent to zero's trick
  }
  
  #Survival (phi must have M x n.primary - 1 dimension for custom updates to work)
  #without individual or primary occasion effects, use for loop to plug into phi[i,g]
  phi.fixed ~ dunif(0,1)
  for(i in 1:M){
    for(g in 1:(n.primary-1)){ #plugging same individual phi's into each primary occasion for custom update
      phi[i,g] <- phi.fixed
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.primary] ~ dSurvival(phi=phi[i,1:(n.primary-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
    #telemetry survival likelihood
    #fixes z states, 1=known alive, 0=known dead, NA=unknown
    #currently assume censoring is uninformative
    tel.z.states[i,1:n.primary] ~ dSurvivalTel(z=z[i,1:n.primary],z.super=z.super[i])
  }
  
  ##Observation Model##
  #sample type observation model priors (Dirichlet), fixed across primary occasions
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])
  sigma.fixed ~ dunif(0,10)
  theta.d.fixed ~ dunif(0.001,500)
  for(g in 1:n.primary){ #sigma informed by data except in primary occasions with no capture effort and no telemetry
    # sigma[g] ~ dunif(0,10) #sigma varies by occasion, shared across methods
    sigma[g] <- sigma.fixed #sigma fixed across primary occasions, shared across methods
  }
  #marking process
  for(g in 1:n.mark.g){
    p0[g] ~ dunif(0,1) #p0 varies by primary occasion
    for(i in 1:M){
      pd[i,mark.g[g],1:J.mark[mark.g[g]]] <- GetDetectionProb(
           s=s[i,1:2],X=X.mark[mark.g[g],1:J.mark[mark.g[g]],1:2],
            J=J.mark[mark.g[g]],sigma=sigma[mark.g[g]],p0=p0[g],z=z[i,mark.g[g]],z.super=z.super[i])
      y.mark[i,mark.g[g],1:J.mark[mark.g[g]]] ~ dBinomialVector(
        pd[i,mark.g[g],1:J.mark[mark.g[g]]],K1D=K1D.mark[mark.g[g],1:J.mark[mark.g[g]]],
        z=z[i,mark.g[g]],z.super=z.super[i])
    }
  }
  #sighting process
  for(g in 1:n.sight.g){
    lam0[g] ~ dunif(0,5) #lam0 varies by primary occasion
    theta.d[g] <- theta.d.fixed
    for(i in 1:M){
      lam[i,sight.g[g],
          1:J.sight[sight.g[g]]] <- GetDetectionRate(
            s=s[i,1:2],X=X.sight[sight.g[g],1:J.sight[sight.g[g]],1:2],
             J=J.sight[sight.g[g]],sigma=sigma[sight.g[g]],lam0=lam0[g],
            z=z[i,sight.g[g]],z.super=z.super[i])
      
    }
    for(i in 1:n.marked.all){ #marked and identified detections
      for(k in 1:K.sight[sight.g[g]]){
        y.mID[i,sight.g[g],1:J.sight[sight.g[g]],k] ~ dPoissonVector(
          lam=lam[i,sight.g[g],1:J.sight[sight.g[g]]]*K2D.sight[sight.g[g],1:J.sight[sight.g[g]],k]*
            theta.marked[1],z=z[i,sight.g[g]],z.super=z.super[i],mark.states=mark.states[i,sight.g[g],k])
      }
    }
    #Unidentified detections by type
    #1 marked with no ID detections
    #sum up lambda contributions of marked individuals when they are marked
    for(k in 1:K.sight[sight.g[g]]){
      bigLam.marked[sight.g[g],1:J.sight[sight.g[g]],k] <- GetbigLam(
        lam=lam[1:n.marked.all,sight.g[g],1:J.sight[sight.g[g]]],
        z=z[1:n.marked.all,sight.g[g]]*mark.states[1:n.marked.all,sight.g[g],k])
      
      lam.mnoID[sight.g[g],1:J.sight[sight.g[g]],k] <- bigLam.marked[sight.g[g],1:J.sight[sight.g[g]],k]*
        K2D.sight[sight.g[g],1:J.sight[sight.g[g]],k]*theta.marked[2]
      
      y.mnoID[sight.g[g],1:J.sight[sight.g[g]],k] ~ dPoissonVector(
        lam.mnoID[sight.g[g],1:J.sight[sight.g[g]],k],z=1,z.super=1,mark.states=1) #plug in z,z.super,mark.states=1 to reuse dPoissonVector
    }
    
    #2 unmarked detections
    #sum up lambda contributions of always unmarked individuals and marked individuals in primary occasions not marked
    for(k in 1:K.sight[sight.g[g]]){
      bigLam.unmarked[sight.g[g],1:J.sight[sight.g[g]],k] <- GetbigLam(
        lam=lam[1:M,sight.g[g],1:J.sight[sight.g[g]]],
        z=z.super[1:M]*z[1:M,sight.g[g]]*(1-mark.states[1:M,sight.g[g],k]))
      
      lam.um[sight.g[g],1:J.sight[sight.g[g]],k] <- bigLam.unmarked[sight.g[g],1:J.sight[sight.g[g]],k]*
        K2D.sight[sight.g[g],1:J.sight[sight.g[g]],k]*theta.unmarked[2]
      
      y.um[sight.g[g],1:J.sight[sight.g[g]],k] ~ dPoissonVector(
        lam.um[sight.g[g],1:J.sight[sight.g[g]],k],
        z=1,z.super=1,mark.states=1) #plug in z,z.super,mark.states=1 to reuse dPoissonVector
    }
    #3 unknown marked status
    for(k in 1:K.sight[sight.g[g]]){
      lam.unk[sight.g[g],1:J.sight[sight.g[g]],k] <- bigLam.marked[sight.g[g],1:J.sight[sight.g[g]],k]*
        K2D.sight[sight.g[g],1:J.sight[sight.g[g]],k]*theta.marked[3] + 
        bigLam.unmarked[sight.g[g],1:J.sight[sight.g[g]],k]*
        K2D.sight[sight.g[g],1:J.sight[sight.g[g]],k]*theta.unmarked[3]
      
      y.unk[sight.g[g],1:J.sight[sight.g[g]],k] ~ dPoissonVector(
        lam.unk[sight.g[g],1:J.sight[sight.g[g]],k],
        z=1,z.super=1,mark.states=1) #plug in z,z.super,mark.states=1 to reuse dPoissonVector
    }
    #integrated detector-session gamma random effects
    trap.RE.dummy[sight.g[g]] ~ dDetectorNBCorrection(
      y.trap.total=y.trap.total[sight.g[g],1:J.sight[sight.g[g]],1:K.sight[sight.g[g]]],
      bigLam.marked=bigLam.marked[sight.g[g],1:J.sight[sight.g[g]],1:K.sight[sight.g[g]]],
      bigLam.unmarked=bigLam.unmarked[sight.g[g],1:J.sight[sight.g[g]],1:K.sight[sight.g[g]]],
      K2D=K2D.sight[sight.g[g],1:J.sight[sight.g[g]],1:K.sight[sight.g[g]]],
      theta=theta.d[g],J=J.sight[sight.g[g]],K=K.sight[sight.g[g]])
  }
  #Telemetry informs activity centers and sigma
  for(i in 1:n.tel.inds){
    for(g in 1:n.tel.sessions[i]){
      locs[i,g,1:max.n.tel.locs,1:2] ~ dNormVector(s=s[tel.ID[i],1:2],
                                                   sigma=sigma[tel.session[i,g]],
                                                   n.locs.ind=n.locs.ind[i,g],
                                                   max.n.tel.locs=max.n.tel.locs)
    }
  }
})

#custom updates:
#1) for marked individuals: update z.start, then update z.stop
#2) for unmarked individuals: update entire z vectors
#3) N.super/z.super update with Herliansyah et al. efficiency improvement
#4) s update with Herliansyah et al. efficiency improvement
