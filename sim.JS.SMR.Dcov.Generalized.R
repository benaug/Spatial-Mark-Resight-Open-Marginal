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

sim.JS.SMR.Dcov.Generalized <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
                            phi=NA,gamma=NA,n.year=NA,
                            theta.marked=NA,theta.unmarked=NA,
                            K.mark=NA,K.sight=NA,K1D.mark=NA,K1D.sight=NA,
                            p0=NA,lam0=NA,sigma=NA,X.mark=NA,X.sight=NA,buff=buff,xlim=NA,
                            ylim=NA,res=NA,
                            mark.year.pars=mark.year.pars,mark.protocol=mark.protocol,
                            n.tel.locs=n.tel.locs){
  
  J.mark <- J.sight <- rep(NA,n.year)
  for(g in 1:n.year){
    X.mark[[g]] <- as.matrix(X.mark[[g]])
    X.sight[[g]] <- as.matrix(X.sight[[g]])
    J.mark[g] <- nrow(X.mark[[g]])
    J.sight[g] <- nrow(X.sight[[g]])
  }
  
  #trap operation - marking process
  if(!any(is.na(K1D.mark))){
    if(length(K1D.mark)!=n.year)stop("K1D.mark must be a list of length n.year")
    for(g in 1:n.year){
      if(any(K1D.mark[[g]]>K.mark[g])){
        stop("Some entries in K1D.mark[[g]] are greater than K.mark[g].")
      }
      if(length(K1D.mark[[g]])!=J.mark[g]){
        stop("K1D.mark[[g]] vector must be of length J.mark[g].")
      }
    }
  }else{
    print("K1D.mark not provided, assuming trap operation is perfect.")
    K1D.mark <- vector("list",n.year)
    for(g in 1:n.year){
      K1D.mark[[g]] <- rep(K.mark[g],J.mark[g])
    }
  }
  
  #trap operation - sighting process
  if(!any(is.na(K1D.sight))){
    if(length(K1D.sight)!=n.year)stop("K1D.sight must be a list of length n.year")
    for(g in 1:n.year){
      if(any(K1D.sight[[g]]>K.sight[g])){
        stop("Some entries in K1D.sight[[g]] are greater than K.sight[g].")
      }
      if(length(K1D.sight[[g]])!=J.sight[g]){
        stop("K1D.sight[[g]] vector must be of length J.sight[g].")
      }
    }
    print("K1D.sight not provided, assuming trap operation is perfect.")
    K1D.sight <- rep(K.sight,J.sight)
  }else{
    print("K1D.sight not provided, assuming trap operation is perfect.")
    K1D.sight <- vector("list",n.year)
    for(g in 1:n.year){
      K1D.sight[[g]] <- rep(K.sight[g],J.sight[g])
    }
  }
  
  #Population dynamics
  N <- rep(NA,n.year)
  N.recruit <- N.survive <- ER <- rep(NA,n.year-1)
  #get expected N in year 1 from D.cov parameters
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
  z <- matrix(0,N[1],n.year)
  z[1:N[1],1] <- 1
  for(g in 2:n.year){
    #Simulate recruits
    ER[g-1] <- N[g-1]*gamma[g-1]
    N.recruit[g-1] <- rpois(1,ER[g-1])
    if(N.recruit[g-1]>0){
      #add recruits to z
      z.dim.old <- nrow(z)
      z <- rbind(z,matrix(0,nrow=N.recruit[g-1],ncol=n.year))
      z[(z.dim.old+1):(z.dim.old+N.recruit[g-1]),g] <- 1
    }
    #Simulate survival
    idx <- which(z[,g-1]==1)
    z[idx,g] <- rbinom(length(idx),1,phi[g-1])
    N.survive[g-1] <- sum(z[,g-1]==1&z[,g]==1)
    N[g] <- N.recruit[g-1]+N.survive[g-1]
  }

  if(any(N.recruit+N.survive!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")

  z.start <- apply(z,1,function(x){which(x==1)[1]})
  z.stop <- n.year-apply(z,1,function(x){which(rev(x)==1)[1]})+1
  
  
  #detection
  J.mark.max <- max(J.mark)
  K.mark.max <- max(K.mark)
  J.sight.max <- max(J.sight)
  K.sight.max <- max(K.sight)

  #simulate activity centers - fixed through time
  N.super <- nrow(z)
  # simulate a population of activity centers
  pi.cell <- lambda.cell/sum(lambda.cell)
  s.cell <- sample(1:n.cells,N.super,prob=pi.cell,replace=TRUE)
  #distribute activity centers uniformly inside cells
  s <- matrix(NA,nrow=N.super,ncol=2)
  for(i in 1:N.super){
    s.xlim <- dSS[s.cell[i],1] + c(-res,res)/2
    s.ylim <- dSS[s.cell[i],2] + c(-res,res)/2
    s[i,1] <- runif(1,s.xlim[1],s.xlim[2])
    s[i,2] <- runif(1,s.ylim[1],s.ylim[2])
  }
  
  #Capture and mark individuals
  y.mark <- pd <- array(0,dim=c(N.super,n.year,J.mark.max))
  for(g in 1:n.year){
    D.mark <- e2dist(s,X.mark[[g]])
    pd[,g,1:J.mark[g]] <- p0[g]*exp(-D.mark*D.mark/(2*sigma[g]*sigma[g]))
    for(i in 1:N.super){
      if(z[i,g]==1){
        y.mark[i,g,1:J.mark[g]] <- rbinom(J.mark[g],size=K1D.mark[[g]],prob=pd[i,g,1:J.mark[g]])
      }
    }
  }
  
  #resight individuals
  lamd <- y <- array(0,dim=c(N.super,n.year,J.sight.max))
  for(g in 1:n.year){
    D <- e2dist(s,X.sight[[g]])
    lamd[,g,1:J.sight[g]] <- lam0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
    for(i in 1:N.super){
      if(z[i,g]==1){
          y[i,g,1:J.sight[g]] <- rpois(J.sight[g],K1D.sight[[g]]*lamd[i,g,1:J.sight[g]])
      }
    }
  }
  
  if(sum(y.mark)==0)stop("No individuals captured. Reconsider parameter settings.")
  if(sum(y)==0)stop("No individuals resighted. Reconsider parameter settings.")

  #store true data for debugging
  truth <- list(y.mark=y.mark,y=y,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,s=s)

  #mark/telemetry data
  #deploy collars to individuals captured in marking process
  mark.caps <- 1*apply(y.mark,c(1,2),sum)
  mark.states <- z*0 #0: unmarked, 1: marked
  #observed data, not true states (because we don't know if dead)
  eligible.states <- matrix(1,N.super,n.year) #eligible based on mark.states collaring history, may be dead and eligible
  for(g in 1:n.year){
    mark.g <- which(mark.caps[,g]>0&eligible.states[,g]==1)
    if(length(mark.g)>0){
      for(i in mark.g){
        mark.life <- rtruncpois(1,lambda=mark.year.pars[1],lower=mark.year.pars[2],upper=mark.year.pars[3])
        end.year <- min(g+mark.life-1,n.year)
        mark.states[i,g:end.year] <- 1
        if(mark.life>1&mark.protocol==1){ #if we don't replace marks on capture, make ineligible
          if(g<n.year){
            eligible.states[i,(g+1):end.year] <- 0
          }
        }
      }
    }
  }
  ID.marked <- vector("list",n.year)
  for(g in 1:n.year){
    ID.marked[[g]] <- which(mark.states[,g]==1)
  }
  
  tel.z.states <- z*NA #0: dead, 1: alive, NA unknown
  tel.z.states[which(mark.states==1&z==1)] <- 1
  tel.z.states[which(mark.states==1&z==0)] <- 0
  #if you observe a death, fill in 0s to the end
  for(i in 1:N.super){
    idx <- which(tel.z.states[i,]==0)
    if(length(idx)>0){
      tel.z.states[i,max(idx):n.year] <- 0
    }
  }
  ID.marked.all <- sort(unique(unlist(ID.marked)))
  n.marked.all <- length(ID.marked.all)
  n.marked <- sapply(ID.marked,length)
  
  #sighting event process V2
  y.event <- array(0,dim=c(N.super,n.year,J.sight.max,3))
  y.mID <- array(0,dim=c(n.marked.all,n.year,J.sight.max))
  y.mnoID <- y.um <- y.unk <- matrix(0,n.year,J.sight.max)
  
  for(g in 1:n.year){
    #loop over cells with positive counts
    idx <- which(y[,g,]>0,arr.ind=TRUE)
    for(l in 1:nrow(idx)){
      if(mark.states[idx[l,1],g]==1){ #if marked
        y.event[idx[l,1],g,idx[l,2],] <- rmultinom(1,y[idx[l,1],g,idx[l,2]],theta.marked)
      }else{#if unmarked
        y.event[idx[l,1],g,idx[l,2],] <- rmultinom(1,y[idx[l,1],g,idx[l,2]],c(0,theta.unmarked,1-theta.unmarked))
      }
    }
    marked.inds <- which(mark.states[,g]==1)
    unmarked.inds <- which(mark.states[,g]==0)
    y.mID[,g,] <- apply(y.event[ID.marked.all,g,,1],c(1,2),sum) #include all marked individuals for consistent individual numbers across years
    if(n.marked[g]>0){
      if(n.marked[g]==1){
        y.mnoID[g,] <- y.event[marked.inds,g,,2]
        y.unk[g,] <- y.event[marked.inds,g,,3] + apply(y.event[unmarked.inds,g,,3],2,sum)
      }else{
        y.mnoID[g,] <- apply(y.event[marked.inds,g,,2],2,sum)
        y.unk[g,] <- apply(y.event[marked.inds,g,,3],2,sum) + apply(y.event[unmarked.inds,g,,3],2,sum)
      }
    }else{
      y.mnoID[g,] <- rep(0,J.sight[g])
      y.unk[g,] <- apply(y.event[unmarked.inds,,3],2,sum) #no marked counts to add
    }
    y.um[g,] <- apply(y.event[unmarked.inds,g,,2],2,sum)
    if(!sum(y[,g,])==(sum(y.mID[,g,])+sum(y.mnoID[g,])+sum(y.um[g,])+sum(y.unk[g,])))stop("data simulator bug")
  }
  
  #old sighting event process
  # y.mID <- y.mnoID <- y.um <- y.unk <- vector("list",n.year)
  # for(g in 1:n.year){
  #   y.event <- array(0,dim=c(N.super,J.sight[g],3))
  #   #loop over cells with positive counts
  #   idx <- which(y[,g,]>0,arr.ind=TRUE)
  #   for(l in 1:nrow(idx)){
  #     if(mark.states[idx[l,1],g]==1){ #if marked
  #       y.event[idx[l,1],idx[l,2],] <- rmultinom(1,y[idx[l,1],g,idx[l,2]],theta.marked)
  #     }else{#if unmarked
  #       y.event[idx[l,1],idx[l,2],] <- rmultinom(1,y[idx[l,1],g,idx[l,2]],c(0,theta.unmarked,1-theta.unmarked))
  #     }
  #   }
  #   marked.inds <- which(mark.states[,g]==1)
  #   unmarked.inds <- which(mark.states[,g]==0)
  #   y.mID[[g]] <- apply(y.event[ID.marked.all,,1],c(1,2),sum) #include all marked individuals for consistent individual numbers across years
  #   if(n.marked[g]>0){
  #     if(n.marked[g]==1){
  #       y.mnoID[[g]] <- y.event[marked.inds,,2]
  #       y.unk[[g]] <- y.event[marked.inds,,3] + apply(y.event[unmarked.inds,,3],2,sum)
  #     }else{
  #       y.mnoID[[g]] <- apply(y.event[marked.inds,,2],2,sum)
  #       y.unk[[g]] <- apply(y.event[marked.inds,,3],2,sum) + apply(y.event[unmarked.inds,,3],2,sum)
  #     }
  #   }else{
  #     y.mnoID[[g]] <- rep(0,J.sight[g])
  #     y.unk[[g]] <- apply(y.event[unmarked.inds,,3],2,sum) #no marked counts to add
  #   }
  #   y.um[[g]] <- apply(y.event[unmarked.inds,,2],2,sum)
  #   if(!sum(y[,g,])==(sum(y.mID[[g]])+sum(y.mnoID[[g]])+sum(y.um[[g]])+sum(y.unk[[g]])))stop("data simulator bug")
  # }
  
  mark.states <- mark.states[ID.marked.all,]
  tel.z.states <- tel.z.states[ID.marked.all,]
  #renumber ID.marked and ID.marked.all in new order after discarding unmarked guys in numbering
  #reorder y, z, s first
  ID.unmarked.all <- setdiff(1:N.super,ID.marked.all)
  s <- s[c(ID.marked.all,ID.unmarked.all),]
  z <- z[c(ID.marked.all,ID.unmarked.all),]
  y.mark <- y.mark[c(ID.marked.all,ID.unmarked.all),,]
  y <- y[c(ID.marked.all,ID.unmarked.all),,]
  #update truth
  truth$s <- s
  truth$z <- z
  truth$y <- y
  
  #reorder marked guys
  for(g in 1:n.year){
    ID.marked[[g]] <- which(mark.states[,g]==1)
  }
  ID.marked.all <- 1:n.marked.all
  
  #discard uncaptured individuals in marking process. keep all marked individuals across years
  # y.mark.list <- vector("list",n.year)
  # for(g in 1:n.year){
  #   y.mark.list[[g]] <- y.mark[ID.marked.all,g,1:J.mark[g]]
  # }
  y.mark <- y.mark[ID.marked.all,,]
  
  
  #Telemetry observations for marked individuals. using all marked individuals here
  tel.inds <- ID.marked
  if(n.tel.locs>0){
    locs <- array(NA,dim=c(max(n.marked),n.year,max(n.tel.locs),2))
    for(g in 1:n.year){
      if(n.marked[g]>0){
        for(i in 1:n.marked[g]){
          for(j in 1:n.tel.locs){
            locs[i,g,j,] <- c(rnorm(1,s[tel.inds[[g]][i],1],sigma),rnorm(1,s[tel.inds[[g]][i],2],sigma))
          }
        }
      }
    }
  }else{
    locs <- NA
  }
  
  return(list(y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk, #observed data
              n.year=n.year,n.marked=n.marked,n.marked.all=n.marked.all,locs=locs,
              ID.marked=ID.marked,ID.marked.all=ID.marked.all,
              mark.states=mark.states,tel.z.states=tel.z.states,
              N=N,N.recruit=N.recruit,N.survive=N.survive,N.super=N.super,X.mark=X.mark,X.sight=X.sight,
              K.mark=K.mark,K.sight=K.sight,
              J.mark=J.mark,J.sight=J.sight,K1D.mark=K1D.mark,K1D.sight=K1D.sight,tel.inds=tel.inds,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.y1=lambda.y1,
              truth=truth))
}
