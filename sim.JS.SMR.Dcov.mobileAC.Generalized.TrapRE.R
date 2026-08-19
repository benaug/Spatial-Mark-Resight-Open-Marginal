e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

rtruncpois <- function(n,lambda,lower=0,upper=Inf){
  p.lo <- ppois(lower-1,lambda)
  p.hi <- ppois(upper,lambda)
  u <- runif(n,min=p.lo,max=p.hi)
  qpois(u,lambda)
}

sim.JS.SMR.Dcov.mobileAC.Generalized.TrapRE <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
                                                        phi=NA,gamma=NA,n.primary=NA,
                                                        theta.marked=NA,theta.unmarked=NA,
                                                        K.mark=NA,K.sight=NA,K1D.mark=NA,K1D.sight=NA,
                                                        p0=NA,lam0=NA,theta.d=NA,sigma=NA,sigma.move=NA,rsf.beta=NA,
                                                        X.mark=NA,X.sight=NA,buff=buff,xlim=NA,
                                                        ylim=NA,res=NA,
                                                        mark.g.pars=NA,mark.protocol=NA,
                                                        n.tel.locs=NA,p.mark=NA){
  
  J.mark <- J.sight <- rep(NA,n.primary)
  for(g in 1:n.primary){
    X.mark[[g]] <- as.matrix(X.mark[[g]])
    X.sight[[g]] <- as.matrix(X.sight[[g]])
    J.mark[g] <- nrow(X.mark[[g]])
    J.sight[g] <- nrow(X.sight[[g]])
    if((K.mark[g]==0)!=(J.mark[g]==0)){
      stop(paste("X.mark and K.mark inconsistent in session",g))
    }
    if((K.sight[g]==0)!=(J.sight[g]==0)){
      stop(paste("X.sight and K.sight inconsistent in session",g))
    }
  }
  
  #trap operation - marking process
  if(!any(is.na(K1D.mark))){
    if(length(K1D.mark)!=n.primary)stop("K1D.mark must be a list of length n.primary")
    for(g in 1:n.primary){
      if(any(K1D.mark[[g]]>K.mark[g])){
        stop("Some entries in K1D.mark[[g]] are greater than K.mark[g].")
      }
      if(length(K1D.mark[[g]])!=J.mark[g]){
        stop("K1D.mark[[g]] vector must be of length J.mark[g].")
      }
    }
  }else{
    print("K1D.mark not provided, assuming trap operation is perfect.")
    K1D.mark <- vector("list",n.primary)
    for(g in 1:n.primary){
      K1D.mark[[g]] <- rep(K.mark[g],J.mark[g])
    }
  }
  
  #trap operation - sighting process
  if(!any(is.na(K1D.sight))){
    if(length(K1D.sight)!=n.primary)stop("K1D.sight must be a list of length n.primary")
    for(g in 1:n.primary){
      if(any(K1D.sight[[g]]>K.sight[g])){
        stop("Some entries in K1D.sight[[g]] are greater than K.sight[g].")
      }
      if(length(K1D.sight[[g]])!=J.sight[g]){
        stop("K1D.sight[[g]] vector must be of length J.sight[g].")
      }
    }
  }else{
    print("K1D.sight not provided, assuming trap operation is perfect.")
    K1D.sight <- vector("list",n.primary)
    for(g in 1:n.primary){
      K1D.sight[[g]] <- rep(K.sight[g],J.sight[g])
    }
  }
  if(length(theta.d)!=n.primary){
    stop("theta.d must have length n.primary.")
  }
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      if(is.na(theta.d[g])|theta.d[g]<=0){
        stop("theta.d must be positive for every sighting session.")
      }
    }
  }
  
  #Population dynamics
  N <- rep(NA,n.primary)
  N.recruit <- N.survive <- ER <- rep(NA,n.primary-1)
  #get expected N in primary occasion 1 from D.cov parameters
  cellArea <- res^2
  lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
  lambda.y1 <- sum(lambda.cell)
  N[1] <- rpois(1,lambda.y1)
  
  #recreate some Dcov things so we can pass fewer arguments into this function
  x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
  y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
  dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
  cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
  n.cells <- nrow(dSS)
  n.cells.x <- length(x.vals)
  n.cells.y <- length(y.vals)
  
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z <- matrix(0,N[1],n.primary)
  z[1:N[1],1] <- 1
  for(g in 2:n.primary){
    #Simulate recruits
    ER[g-1] <- N[g-1]*gamma[g-1]
    N.recruit[g-1] <- rpois(1,ER[g-1])
    if(N.recruit[g-1]>0){
      #add recruits to z
      z.dim.old <- nrow(z)
      z <- rbind(z,matrix(0,nrow=N.recruit[g-1],ncol=n.primary))
      z[(z.dim.old+1):(z.dim.old+N.recruit[g-1]),g] <- 1
    }
    #Simulate survival
    idx <- which(z[,g-1]==1)
    z[idx,g] <- rbinom(length(idx),1,phi[g-1])
    N.survive[g-1] <- sum(z[,g-1]==1&z[,g]==1)
    N[g] <- N.recruit[g-1]+N.survive[g-1]
  }
  
  if(any(N.recruit+N.survive!=N[2:n.primary]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  z.start <- apply(z,1,function(x){which(x==1)[1]})
  z.stop <- n.primary-apply(z,1,function(x){which(rev(x)==1)[1]})+1
  
  
  #detection
  J.mark.max <- max(J.mark)
  K.mark.max <- max(K.mark)
  J.sight.max <- max(J.sight)
  K.sight.max <- max(K.sight)
  
  # simulate a population of activity centers for primary occasion 1 proportional to D.cov
  N.super <- nrow(z)
  library(truncnorm)
  pi.cell <- lambda.cell/sum(lambda.cell)
  s.cell <- matrix(NA,N.super,n.primary)
  s.cell[,1] <- sample(1:n.cells,N.super,prob=pi.cell,replace=TRUE)
  #distribute activity centers uniformly inside cells
  s <- array(NA,dim=c(N.super,n.primary,2))
  for(i in 1:N.super){
    s.xlim <- dSS[s.cell[i,1],1] + c(-res,res)/2
    s.ylim <- dSS[s.cell[i,1],2] + c(-res,res)/2
    s[i,1,1] <- runif(1,s.xlim[1],s.xlim[2])
    s[i,1,2] <- runif(1,s.ylim[1],s.ylim[2])
  }
  #subsequent primary occasions
  avail.dist <- use.dist <- array(NA,dim=c(N.super,n.primary-1,n.cells))
  rsf <- exp(rsf.beta*D.cov)
  rsf[InSS==0] <- 0 #disallow individuals moving into nonhabitat
  for(g in 2:n.primary){
    for(i in 1:N.super){
      avail.dist[i,g-1,] <- getAvail(s=s[i,g-1,1:2],sigma=sigma.move,res=res,x.vals=x.vals,
                                     y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
      use.dist[i,g-1,] <- rsf*avail.dist[i,g-1,]
      use.dist[i,g-1,] <- use.dist[i,g-1,]/sum(use.dist[i,g-1,])
      #move AC - select new cell
      s.cell[i,g] <- sample(1:n.cells,1,replace=TRUE,prob=use.dist[i,g-1,])
      #choose location inside cell
      s.xlim <- dSS[s.cell[i,g],1] + c(-res,res)/2
      s.ylim <- dSS[s.cell[i,g],2] + c(-res,res)/2
      #choose new location inside cell
      s[i,g,1] <- rtruncnorm(1,a=s.xlim[1],b=s.xlim[2],mean=s[i,g-1,1],sd=sigma.move)
      s[i,g,2] <- rtruncnorm(1,a=s.ylim[1],b=s.ylim[2],mean=s[i,g-1,2],sd=sigma.move)
    }
  }
  
  #Capture and mark individuals
  y.mark <- pd <- array(0,dim=c(N.super,n.primary,J.mark.max))
  for(g in 1:n.primary){
    if(J.mark[g]>0){
      D.mark <- e2dist(s[,g,],X.mark[[g]])
      pd[,g,1:J.mark[g]] <- p0[g]*exp(-D.mark*D.mark/(2*sigma[g]*sigma[g]))
      for(i in 1:N.super){
        if(z[i,g]==1){
          y.mark[i,g,1:J.mark[g]] <- rbinom(J.mark[g],size=K1D.mark[[g]],prob=pd[i,g,1:J.mark[g]])
        }
      }
    }
  }
  
  #resight individuals with shared detector-session random effects
  lamd <- y <- array(0,dim=c(N.super,n.primary,J.sight.max))
  trap.effect <- matrix(1,nrow=n.primary,ncol=J.sight.max)
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      D <- e2dist(s[,g,],X.sight[[g]])
      lamd[,g,1:J.sight[g]] <- lam0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
      #shared effect across all individuals at detector j in session g
      trap.effect[g,1:J.sight[g]] <- rgamma(J.sight[g],shape=theta.d[g],rate=theta.d[g])
      for(i in 1:N.super){
        if(z[i,g]==1){
          y[i,g,1:J.sight[g]] <- rpois(J.sight[g],K1D.sight[[g]]*lamd[i,g,1:J.sight[g]]*trap.effect[g,1:J.sight[g]])
        }
      }
    }
  }
  
  if(sum(y)==0)stop("No individuals resighted. Reconsider parameter settings.")
  
  #expected proportion of realized N in cell 
  pi.cell <- array(NA,dim=c(n.primary,n.cells))
  pi.cell[1,] <- lambda.cell/sum(lambda.cell)
  for(g in 2:n.primary){
    pi.cell[g,] <- apply(use.dist[z[,g]==1,g-1,,drop=FALSE],3,sum)*InSS
    pi.cell[g,] <- pi.cell[g,]/sum(pi.cell[g,])
  }
  
  #store true data for debugging
  truth <- list(y.mark=y.mark,y=y,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,N.super=N.super,
                s=s,s.cell=s.cell,trap.effect=trap.effect,
                pi.cell=pi.cell,avail.dist=avail.dist,use.dist=use.dist)
  
  #mark/telemetry data
  #deploy collars to individuals captured in marking process
  mark.caps <- 1*apply(y.mark,c(1,2),sum)
  ID.cap.all <- sort(unique(which(rowSums(mark.caps)>0)))
  n.cap.all <- length(ID.cap.all)
  mark.deploy <- matrix(0,N.super,n.primary) #actual marks deployed only
  mark.states <- z*0 #0: unmarked, 1: marked
  tel.z.states <- z*NA
  #observed data, not true states (because we don't know if dead)
  eligible.states <- matrix(1,N.super,n.primary) #eligible based on mark.states collaring history, may be dead and eligible
  for(g in 1:n.primary){
    cap.g <- which(mark.caps[,g]>0&eligible.states[,g]==1)
    if(length(cap.g)>0){
      deploy.g <- rbinom(length(cap.g),1,p.mark[g])
      mark.deploy[cap.g,g] <- deploy.g
      mark.g <- cap.g[which(deploy.g==1)]
      if(length(mark.g)>0){
        for(i in mark.g){
          mark.life <- rtruncpois(1,lambda=mark.g.pars[1],lower=mark.g.pars[2],upper=mark.g.pars[3])
          end.g <- min(g+mark.life-1,n.primary)
          mark.states[i,g:end.g] <- 1
          tel.z.states[i,g:end.g] <- 1
          if(mark.life>1&mark.protocol==1){ #if we don't replace marks on capture, make ineligible
            if(g<n.primary){
              eligible.states[i,(g+1):end.g] <- 0
            }
          }
        }
      }
    }
  }
  if(sum(mark.deploy)==0)stop("No individuals marked. Reconsider parameter settings.")
  
  #switch states to observed deaths when z==0 
  tel.z.states[which(tel.z.states==1&z==0)] <- 0
  mark.states[which(mark.states==1&z==0)] <- 0
  ID.marked <- vector("list",n.primary)
  for(g in 1:n.primary){
    ID.marked[[g]] <- which(mark.states[,g]==1)
  }
  #if you observe a death, fill in 0s to the end
  for(i in 1:N.super){
    idx <- which(tel.z.states[i,]==0)
    if(length(idx)>0){
      tel.z.states[i,max(idx):n.primary] <- 0
    }
  }
  ID.marked.all <- sort(unique(unlist(ID.marked)))
  n.marked.all <- length(ID.marked.all)
  n.marked <- sapply(ID.marked,length)
  
  #sighting event process V2
  y.event <- array(0,dim=c(N.super,n.primary,J.sight.max,3))
  y.mID <- array(0,dim=c(n.marked.all,n.primary,J.sight.max))
  y.mnoID <- y.um <- y.unk <- matrix(0,n.primary,J.sight.max)
  
  for(g in 1:n.primary){
    if(K.sight[g]>0){ #skip if no effort in this primary occasion
      #loop over cells with positive counts
      y.g <- matrix(y[,g,1:J.sight[g]],nrow=N.super,ncol=J.sight[g])
      idx <- which(y.g>0,arr.ind=TRUE)
      if(nrow(idx)>0){
        for(l in 1:nrow(idx)){
          if(mark.states[idx[l,1],g]==1){ #if marked
            y.event[idx[l,1],g,idx[l,2],] <- rmultinom(1,y[idx[l,1],g,idx[l,2]],theta.marked)
          }else{#if unmarked
            y.event[idx[l,1],g,idx[l,2],] <- rmultinom(1,y[idx[l,1],g,idx[l,2]],c(0,theta.unmarked,1-theta.unmarked))
          }
        }
      }
      marked.inds <- which(mark.states[,g]==1)
      unmarked.inds <- which(mark.states[,g]==0)
      for(i2 in 1:n.marked.all){
        y.mID[i2,g,1:J.sight[g]] <- y.event[ID.marked.all[i2],g,1:J.sight[g],1]
      }
      for(j in 1:J.sight[g]){
        if(length(marked.inds)>0){
          y.mnoID[g,j] <- sum(y.event[marked.inds,g,j,2])
          y.unk[g,j] <- sum(y.event[marked.inds,g,j,3])
        }
        if(length(unmarked.inds)>0){
          y.um[g,j] <- sum(y.event[unmarked.inds,g,j,2])
          y.unk[g,j] <- y.unk[g,j]+sum(y.event[unmarked.inds,g,j,3])
        }
      }
      if(!sum(y[,g,])==(sum(y.mID[,g,])+sum(y.mnoID[g,])+sum(y.um[g,])+sum(y.unk[g,])))stop("data simulator bug")
    }
  }
  
  #simulate telemetry locations for all collared primary occasions
  if(n.tel.locs>0&sum(mark.states)>0){
    n.tel.sessions.vec <- rowSums(tel.z.states==1,na.rm=TRUE)
    tel.ID <- which(n.tel.sessions.vec>0)
    n.tel.inds <- length(tel.ID)
    if(n.tel.inds>0){
      n.tel.sessions <- n.tel.sessions.vec[tel.ID] #length n.tel.inds
      max.n.tel.sessions <- max(n.tel.sessions)
      locs <- array(NA,dim=c(n.tel.inds,max.n.tel.sessions,n.tel.locs,2))
      n.locs.ind <- matrix(0,n.tel.inds,max.n.tel.sessions)
      tel.session <- matrix(NA,n.tel.inds,max.n.tel.sessions)
      for(i in 1:n.tel.inds){
        collared.gs <- which(tel.z.states[tel.ID[i],]==1)
        tel.session[i,1:length(collared.gs)] <- collared.gs
        for(gy in 1:length(collared.gs)){
          g <- collared.gs[gy]
          locs[i,gy,1:n.tel.locs,1] <- rnorm(n.tel.locs,s[tel.ID[i],g,1],sigma[g])
          locs[i,gy,1:n.tel.locs,2] <- rnorm(n.tel.locs,s[tel.ID[i],g,2],sigma[g])
          n.locs.ind[i,gy] <- n.tel.locs
        }
      }
    }else{
      locs <- tel.ID <- tel.session <- NA
      n.tel.inds <- 0
      n.locs.ind <- NA
      n.tel.sessions <- NA
    }
  }else{
    print("no individuals captured, no telemetry")
    locs <- tel.ID <- tel.session <- NA
    n.tel.inds <- 0
    n.locs.ind <- NA
    n.tel.sessions <- NA
  }
  tel.ID.g <- vector("list",n.primary)
  for(g in 1:n.primary){
    collared.g <- which(tel.z.states[,g]==1)
    if(length(collared.g)>0){
      tel.ID.g[[g]] <- collared.g
    }
  }
  
  #simulate telemetry locations
  # if(n.tel.locs>0&sum(y.mark)>0){
  #   n.tel.sessions <- rowSums(tel.z.states==1,na.rm=TRUE)
  #   n.tel.sessions <- n.tel.sessions[ID.marked.all]
  #   n.tel.inds <- sum(n.tel.sessions>0)
  #   tel.session <- matrix(NA,n.tel.inds,n.primary)
  #   max.n.tel.sessions <- max(n.tel.sessions)
  #   locs <- array(NA,dim=c(n.tel.inds,max.n.tel.sessions,n.tel.locs,2))
  #   for(i in 1:n.tel.inds){
  #     tel.session[i,1:n.tel.sessions[i]] <- which(tel.z.states[ID.marked.all[i],]==1)
  #     for(g in 1:n.tel.sessions[i]){
  #       #if adding movement, reference correct s primary occasions
  #       locs[i,g,,] <- c(rnorm(n.tel.locs,s[ID.marked.all[i],1],sigma[tel.session[i,g]]),
  #                        rnorm(n.tel.locs,s[ID.marked.all[i],2],sigma[tel.session[i,g]]))
  #     }
  #   }
  #   n.locs.ind <- apply(!is.na(locs[,,,1]),c(1,2),sum)
  #   if(dim(locs)[2]==1){
  #     n.locs.ind <- matrix(rowSums(n.locs.ind),ncol=1)
  #   }
  # }else{
  #   print("no individuals captured, no telemetry")
  #   locs <- tel.session <- NA
  #   n.tel.inds <- 0
  #   n.tel.sessions <- NA
  #   n.locs.ind <- NA
  # }
  
  mark.states <- mark.states[ID.marked.all,,drop=FALSE]
  tel.z.states <- tel.z.states[ID.marked.all,,drop=FALSE]
  
  #renumber ID.marked and ID.marked.all in new order after discarding unmarked guys in numbering
  #reorder y, z, s first
  # ID.unmarked.all <- setdiff(1:N.super,ID.marked.all)
  # s <- s[c(ID.marked.all,ID.unmarked.all),]
  # z <- z[c(ID.marked.all,ID.unmarked.all),]
  # y.mark <- y.mark[c(ID.marked.all,ID.unmarked.all),,]
  # y <- y[c(ID.marked.all,ID.unmarked.all),,]
  ID.cap.unmarked.all <- setdiff(ID.cap.all,ID.marked.all)
  ID.unobserved.all <- setdiff(1:N.super,c(ID.marked.all,ID.cap.unmarked.all))
  ID.order <- c(ID.marked.all,ID.cap.unmarked.all,ID.unobserved.all)
  s <- s[ID.order,,,drop=FALSE]
  s.cell <- s.cell[ID.order,,drop=FALSE]
  avail.dist <- avail.dist[ID.order,,,drop=FALSE]
  use.dist <- use.dist[ID.order,,,drop=FALSE]
  z <- z[ID.order,,drop=FALSE]
  y.mark <- y.mark[ID.order,,,drop=FALSE]
  y <- y[ID.order,,,drop=FALSE]
  mark.caps <- mark.caps[ID.order,,drop=FALSE]
  mark.deploy <- mark.deploy[ID.order,,drop=FALSE]
  #update truth
  truth$s <- s
  truth$s.cell <- s.cell
  truth$avail.dist <- avail.dist
  truth$use.dist <- use.dist
  truth$z <- z
  truth$y <- y
  truth$y.mark <- y.mark
  
  #reorder marked guys
  for(g in 1:n.primary){
    ID.marked[[g]] <- which(mark.states[,g]==1)
  }
  tel.ID <- match(tel.ID,ID.marked.all)
  for(g in 1:n.primary){
    if(length(tel.ID.g[[g]])>0){
      tel.ID.g[[g]] <- match(tel.ID.g[[g]],ID.marked.all)
    }
  }
  ID.marked.all <- 1:n.marked.all
  
  #discard uncaptured individuals in marking process. keep marked and unmarked (no mark deployed) captured in marking process
  y.mark <- y.mark[1:n.cap.all,,,drop=FALSE]
  ID.cap.all <- 1:n.cap.all
  mark.caps <- mark.caps[1:n.cap.all,,drop=FALSE]
  mark.deploy <- mark.deploy[1:n.cap.all,,drop=FALSE]
  
  return(list(y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk, #observed data
              n.primary=n.primary,n.marked=n.marked,n.marked.all=n.marked.all,
              ID.cap.all=ID.cap.all,n.cap.all=n.cap.all,mark.deploy=mark.deploy,mark.caps=mark.caps,
              locs=locs,n.tel.inds=n.tel.inds,n.tel.sessions=n.tel.sessions,tel.session=tel.session,
              n.locs.ind=n.locs.ind,tel.ID=tel.ID,tel.ID.g=tel.ID.g,
              ID.marked=ID.marked,ID.marked.all=ID.marked.all,
              mark.states=mark.states,tel.z.states=tel.z.states,
              N=N,N.recruit=N.recruit,N.survive=N.survive,N.super=N.super,X.mark=X.mark,X.sight=X.sight,
              K.mark=K.mark,K.sight=K.sight,
              J.mark=J.mark,J.sight=J.sight,K1D.mark=K1D.mark,K1D.sight=K1D.sight,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,s=s,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,lambda.y1=lambda.y1,
              truth=truth))
}