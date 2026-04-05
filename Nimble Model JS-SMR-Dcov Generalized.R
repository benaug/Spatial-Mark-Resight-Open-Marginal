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
  for(g in 2:n.year){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly abundance
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.year-1)]) #size of superpopulation
  
  #Recruitment
  gamma.fixed ~ dunif(0,2)
  for(g in 1:(n.year-1)){
    # gamma[g] ~ dunif(0,2) # yearly recruitment priors
    gamma[g] <- gamma.fixed
    ER[g] <- N[g]*gamma[g] #yearly expected recruits
    N.recruit[g] ~ dpois(ER[g]) #yearly realized recruits
  }
  
  #Individual covariates
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]]) #categorical likelihood for this cell, equivalent to zero's trick
  }
  
  #Survival (phi must have M x n.year - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  phi.fixed ~ dunif(0,1)
  for(i in 1:M){
    for(g in 1:(n.year-1)){ #plugging same individual phi's into each year for custom update
      phi[i,g] <- phi.fixed
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.year] ~ dSurvival(phi=phi[i,1:(n.year-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
    #telemetry survival likelihood
    #fixes z states, 1=known alive, 0=known dead, NA=unknown
    #currently assume censoring is uninformative
    tel.z.states[i,1:n.year] ~ dSurvivalTel(z=z[i,1:n.year],z.super=z.super[i])
  }
  
  ##Detection##
  for(g in 1:n.year){
    p0[g] ~ dunif(0,1) #p0 varies by year
    lam0[g] ~ dunif(0,1) #lam0 varies by year
    sigma[g] ~ dunif(0,10) #sigma varies by year
    #sample type observation model priors (Dirichlet), vary by year
    alpha.marked[g,1] <- 1
    alpha.marked[g,2] <- 1
    alpha.marked[g,3] <- 1
    alpha.unmarked[g,1] <- 1
    alpha.unmarked[g,2] <- 1
    theta.marked[g,1:3] ~ ddirch(alpha.marked[g,1:3])
    theta.unmarked[g,1] <- 0
    theta.unmarked[g,2:3] ~ ddirch(alpha.unmarked[g,1:2])
    for(i in 1:M){
      lam[i,g,1:J.sight[g]] <- GetDetectionRate(s = s[i,1:2], X=X.sight[g,1:J.sight[g],1:2], J=J.sight[g],sigma=sigma[g], lam0=lam0[g],
                                          z=z[i,g],z.super=z.super[i])
      pd[i,g,1:J.mark[g]] <- GetDetectionProb(s=s[i,1:2],X=X.mark[g,1:J.mark[g],1:2],J=J.mark[g],sigma=sigma[g],p0=p0[g],
                                              z=z[i,g],z.super=z.super[i])
      y.mark[i,g,1:J.mark[g]] ~ dBinomialVector(pd[i,g,1:J.mark[g]],K1D=K1D.mark[g,1:J.mark[g]],z=z[i,g],z.super=z.super[i])
    }
    for(i in 1:n.marked.all){
      y.mID[i,g,1:J.sight[g]] ~ dPoissonVector(lam=lam[i,g,1:J.sight[g]]*K1D.sight[g,1:J.sight[g]]*theta.marked[g,1],
                                         z=z[i,g],z.super=z.super[i],mark.states=mark.states[i,g]) #marked and identified detections
    }
    #Unidentified detections by type
    #1 marked with no ID detections
    #sum up lambda contributions of marked individuals when they are marked
    bigLam.marked[g,1:J.sight[g]] <- GetbigLam(lam=lam[1:n.marked.all,g,1:J.sight[g]],
                                         z=z[1:n.marked.all,g]*mark.states[1:n.marked.all,g])
    lam.mnoID[g,1:J.sight[g]] <- bigLam.marked[g,1:J.sight[g]]*K1D.sight[g,1:J.sight[g]]*theta.marked[g,2]
    y.mnoID[g,1:J.sight[g]] ~ dPoissonVector(lam.mnoID[g,1:J.sight[g]],z=1,z.super=1,mark.states=1) #plug in z,z.super=1 to reuse dPoissonVector

    #2 unmarked detections
    #sum up lambda contributions of always unmarked individuals and marked individuals in years not marked
    bigLam.unmarked[g,1:J.sight[g]] <- GetbigLam(lam=lam[1:M,g,1:J.sight[g]],z=z.super[1:M]*z[1:M,g]*(1-mark.states[1:M,g]))
    lam.um[g,1:J.sight[g]] <- bigLam.unmarked[g,1:J.sight[g]]*K1D.sight[g,1:J.sight[g]]*theta.unmarked[g,2]
    y.um[g,1:J.sight[g]] ~ dPoissonVector(lam.um[g,1:J.sight[g]],z=1,z.super=1,mark.states=1) #plug in z,z.super=1 to reuse dPoissonVector

    #3 unknown marked status
    lam.unk[g,1:J.sight[g]] <- bigLam.marked[g,1:J.sight[g]]*K1D.sight[g,1:J.sight[g]]*theta.marked[g,3] +
      bigLam.unmarked[g,1:J.sight[g]]*K1D.sight[g,1:J.sight[g]]*theta.unmarked[g,3]
    y.unk[g,1:J.sight[g]] ~ dPoissonVector(lam.unk[g,1:J.sight[g]],z=1,z.super=1,mark.states=1) #plug in z,z.super=1 to reuse dPoissonVector
    
    #If you have telemetry in every year
    for(i in 1:n.tel.inds[g]){
      for(m in 1:n.locs.ind[i,g]){
        locs[i,g,m,1] ~ dnorm(s[tel.inds[i,g],1],sd=sigma[g])
        locs[i,g,m,2] ~ dnorm(s[tel.inds[i,g],2],sd=sigma[g])
      }
    }
  }
})

#custom updates:
#1) for marked individuals: update z.start, then update z.stop
#2) for unmarked individuals: update entire z vectors
#3) N.super/z.super update with Herliansyah et al. efficiency improvement
#4) s update with Herliansyah et al. efficiency improvement