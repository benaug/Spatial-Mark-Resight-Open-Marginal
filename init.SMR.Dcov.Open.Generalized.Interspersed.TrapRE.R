e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.Open.Generalized.Interspersed.TrapRE <- function(data,inits=NA,M=NA){
  if(M < (data$n.cap.all)+1) stop("M must be larger than the number of captured individuals plus at least one unmarked individual.")
  library(abind)
  n.marked <- data$n.marked
  n.marked.all <- data$n.marked.all
  n.cap.all <- data$n.cap.all
  n.primary <- data$n.primary
  J.mark <- unlist(lapply(data$X.mark,nrow)) #traps per year
  J.sight <- unlist(lapply(data$X.sight,nrow)) #traps per year
  J.mark.max <- max(J.mark)
  J.sight.max <- max(J.sight)
  K.mark <- data$K.mark
  K.sight <- data$K.sight
  K.mark.max <- max(K.mark)
  K.sight.max <- max(K.sight)
  locs <- data$locs
  mark.states2D <- matrix(0,M,n.primary)
  mark.states2D[1:n.marked.all,] <- data$mark.states2D
  mark.states <- array(0,dim=c(M,n.primary,K.sight.max))
  mark.states[1:n.marked.all,,] <- data$mark.states
  #augment tel.z.states, code NA as 2
  tel.z.states <- matrix(2,M,n.primary)
  tel.z.states[1:n.marked.all,] <- data$tel.z.states
  tel.z.states[is.na(tel.z.states)] <- 2  #use 2 to code NA for nimble function

  #augment y.mark and y.mnoID, pull out y.mnoID, y.um, y.unk
  y.mark <- array(0,dim=c(M,n.primary,J.mark.max))
  y.mark[1:n.cap.all,,] <- data$y.mark
  y.mID <- array(0,dim=c(n.marked.all,n.primary,J.sight.max,K.sight.max))
  y.mID[1:n.marked.all,,,] <- data$y.mID
  y.mnoID <- data$y.mnoID
  y.um <- data$y.um
  y.unk <- data$y.unk
  
  y.trap.total <- array(0,dim=c(n.primary,J.sight.max,K.sight.max))
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      y.trap.total[g,1:J.sight[g],1:K.sight[g]] <-
        apply(y.mID[,g,1:J.sight[g],1:K.sight[g],drop=FALSE],c(3,4),sum)+
        y.mnoID[g,1:J.sight[g],1:K.sight[g]]+
        y.um[g,1:J.sight[g],1:K.sight[g]]+
        y.unk[g,1:J.sight[g],1:K.sight[g]]
    }
  }
  trap.RE.dummy <- rep(0,n.primary)
  
  #reformat these
  ID.marked <- matrix(0,max(n.marked),n.primary)
  X.mark <- array(0,dim=c(n.primary,J.mark.max,2))
  K1D.mark <- matrix(0,n.primary,J.mark.max)
  X.sight <- array(0,dim=c(n.primary,J.sight.max,2))
  K2D.sight <- array(0,dim=c(n.primary,J.sight.max,K.sight.max))
  for(g in 1:n.primary){
    if(n.marked[g]>0){
      ID.marked[1:n.marked[g],g] <- data$ID.marked[[g]]
    }
    if(J.mark[g]>0){
      X.mark[g,1:J.mark[g],1:2] <- data$X.mark[[g]]
      K1D.mark[g,1:J.mark[g]] <- data$K1D.mark[[g]]
    }
    if(J.sight[g]>0){
      X.sight[g,1:J.sight[g],1:2] <- data$X.sight[[g]]
      K2D.sight[g,1:J.sight[g],1:K.sight[g]] <- data$K2D.sight[[g]]
    }
  }
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  p0 <- inits$p0
  lam0 <- inits$lam0
  sigma <- inits$sigma
  theta.d <- inits$theta.d
  
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  
  #but update s.inits for marked individuals and unmarked captured in marking process before assigning latent detections
  has.mark <- apply(y.mark[1:n.cap.all,,],1,sum)>0
  has.mID <- rep(FALSE,n.cap.all)
  has.mID[1:n.marked.all] <- rowSums(y.mID)>0
  has.tel <- (1:n.cap.all)%in%data$tel.ID
  idx <- which(has.mark|has.mID|has.tel)
  
  for(i in idx){
    trps <- matrix(0,nrow=0,ncol=2) #get locations of traps of capture across years for ind i
    for(g in 1:n.primary){
      if(sum(y.mark[i,g,])>0){
        trps.g <- matrix(X.mark[g,which(y.mark[i,g,]>0),,drop=FALSE],ncol=2,byrow=FALSE)
        trps <- rbind(trps,trps.g)
      }
      if(i<=n.marked.all){
        if(sum(y.mID[i,g,,])>0){
          y.mID.ij <- apply(y.mID[i,g,1:J.sight[g],1:K.sight[g],drop=FALSE],3,sum)
          trps.g <- matrix(X.sight[g,which(y.mID.ij>0),,drop=FALSE],ncol=2,byrow=FALSE)
          trps <- rbind(trps,trps.g)
        }
      }
    }
    if(i%in%data$tel.ID){
      these.locs <- matrix(0,nrow=0,ncol=2) #get telemetry locs across years for tel.ID i
      this.tel.ind <- which(data$tel.ID==i)
      for(g in 1:data$n.tel.sessions[this.tel.ind]){
        if(data$n.locs.ind[this.tel.ind,g]>0){
          locs.g <- matrix(data$locs[this.tel.ind,g,1:data$n.locs.ind[this.tel.ind,g],,drop=FALSE],
                           ncol=2,byrow=FALSE)
          these.locs <- rbind(these.locs,locs.g)
        }
      }
      trps <- rbind(trps,these.locs)
    }
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else if(nrow(trps)==1){
      s.init[i,] <- trps
    }
  }
  #in case s.init exactly on state space boundary
  eps <- 1e-5
  s.init[,1][s.init[,1]<xlim[1]] <- xlim[1] + eps
  s.init[,1][s.init[,1]>xlim[2]] <- xlim[2] - eps
  s.init[,2][s.init[,2]<ylim[1]] <- ylim[1] + eps
  s.init[,2][s.init[,2]>ylim[2]] <- ylim[2] - eps
  
  #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
  e2dist  <-  function (x, y){
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  getCell  <-  function(s,res,cells){
    cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
  }
  alldists <- e2dist(s.init,data$dSS)
  alldists[,data$InSS==0] <- Inf
  for(i in 1:M){
    this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
    if(data$InSS[this.cell]==0){
      cands <- alldists[i,]
      new.cell <- which(alldists[i,]==min(alldists[i,]))[1]
      s.init[i,] <- data$dSS[new.cell,]
    }
  }
  
  y.true <- array(0,dim=c(M,n.primary,J.sight.max,K.sight.max))
  y.true[1:n.marked.all,,,] <- y.mID
  
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      D <- e2dist(s.init,X.sight[g,1:J.sight[g],])
      lamd <- lam0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
      for(j in 1:J.sight[g]){
        if(J.sight[g]>0){
          for(k in 1:K.sight[g]){
            marked.inds <- which(mark.states[,g,k]==1)
            unmarked.inds <- which(mark.states[,g,k]==0)
            if(length(marked.inds)>0){
              prob <- lamd[marked.inds,j]*(tel.z.states[marked.inds,g]!=0)
              if(sum(prob)>0){
                prob <- prob/sum(prob)
                y.true[marked.inds,g,j,k] <- y.true[marked.inds,g,j,k]+
                  rmultinom(1,y.mnoID[g,j,k],prob=prob)
              }
            }
            prob <- lamd[,j]*(1-mark.states[,g,k])*(tel.z.states[,g]!=0)
            if(sum(prob)>0){
              prob <- prob/sum(prob)
              y.true[,g,j,k] <- y.true[,g,j,k]+rmultinom(1,y.um[g,j,k],prob=prob)
            }
            prob <- lamd[,j]*(tel.z.states[,g]!=0)
            if(sum(prob)>0){
              prob <- prob/sum(prob)
              y.true[,g,j,k] <- y.true[,g,j,k]+rmultinom(1,y.unk[g,j,k],prob=prob)
            }
          }
        }
      }
    }
  }
  #initialize z, start with observed guys
  z.init <- matrix(0,M,n.primary)
  z.start.init <- z.stop.init <- rep(0,M)
  y.mark2D <- apply(y.mark,c(1,2),sum)
  y.true2D <- 1*((apply(y.true,c(1,2),sum)+y.mark2D)>0)
  y.true2D[tel.z.states==1] <- 1
  z.init <- 1*(y.true2D>0)
  N.super.init <- sum(rowSums(z.init)>0)
  for(i in 1:M){
    det.idx <- which(y.true2D[i,]>0)
    if(length(det.idx)>0){
      z.start.init[i] <- min(det.idx)
      z.stop.init[i] <- max(det.idx)
      z.init[i,z.start.init[i]:z.stop.init[i]] <- 1
    }
  }
  z.super.init <- 1*(z.start.init>0)

  if(any(tel.z.states[z.init==1]==0))stop("At least one z initialized to 1 when tel.z.states=0. Bug in initialization code.")

  #initialize N structures from z.init
  N.init <- colSums(z.init[z.super.init==1,])
  N.survive.init <- N.recruit.init <- rep(NA,n.primary-1)
  for(g in 2:n.primary){
    N.survive.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1)
    N.recruit.init[g-1] <- N.init[g]-N.survive.init[g-1]
  }
  
  #get y2D constraints for z.start and z.stop update
  y.mID2D <- matrix(0,M,n.primary)
  y.mID2D[1:n.marked.all,] <- apply(y.mID,c(1,2),sum)
  y2D <- y.mark2D + y.mID2D
  for(i in 1:n.marked.all){
    idx <- which(tel.z.states[i,]==1)
    if(length(idx)>0){
      y2D[i,idx] <- 1
    }
  }
  
  #check starting logProbs one session at a time
  for(g in 1:n.primary){
    #marking process
    if(J.mark[g]>0){
      D.mark <- e2dist(s.init, X.mark[g,1:J.mark[g],])
      pd <- p0[g]*exp(-D.mark*D.mark/(2*sigma[g]*sigma[g]))
      logProb <- array(0,dim=c(M,J.mark[g]))
      for(i in 1:M){
        for(j in 1:J.mark[g]){
          logProb[i,j] <- dbinom(y.mark[i,g,j],size=K1D.mark[g,j],prob=pd[i,j],log=TRUE)
        }
      }
      if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marking process, year",g))
    }
    #sighting process
    if(J.sight[g]>0){
      D.sight <- e2dist(s.init, X.sight[g,1:J.sight[g],])
      lamd <- lam0[g]*exp(-D.sight*D.sight/(2*sigma[g]*sigma[g]))
      #marked with ID obs
      if(n.marked[g]>0){
        logProb <- array(0,dim=c(n.marked.all,J.sight[g],K.sight[g]))
        for(i in 1:n.marked.all){
          for(j in 1:J.sight[g]){
            for(k in 1:K.sight[g]){
              if(mark.states[i,g,k]==1){
                logProb[i,j,k] <- dpois(y.mID[i,g,j,k],
                                        lamd[i,j]*K2D.sight[g,j,k],
                                        log=TRUE)
              }else if(y.mID[i,g,j,k]>0){
                logProb[i,j,k] <- -Inf
              }
            }
          }
        }
        if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marked with ID observations, year",g))
      }
      for(k in 1:K.sight[g]){
        marked.inds.gk <- which(mark.states[,g,k]==1)
        unmarked.inds.gk <- which(mark.states[,g,k]==0)
        lamd.mnoID <- rep(0,J.sight[g])
        if(length(marked.inds.gk)>0){
          lamd.mnoID <- colSums(lamd[marked.inds.gk,1:J.sight[g],drop=FALSE])
        }
        lamd.um <- rep(0,J.sight[g])
        if(length(unmarked.inds.gk)>0){
          lamd.um <- colSums(lamd[unmarked.inds.gk,1:J.sight[g],drop=FALSE])
        }
        lamd.unk <- colSums(lamd[1:M,1:J.sight[g],drop=FALSE])
        for(j in 1:J.sight[g]){
          if(!is.finite(dpois(y.mnoID[g,j,k],lamd.mnoID[j]*K2D.sight[g,j,k],log=TRUE))){
            stop(paste("Starting observation model likelihood not finite. Marked with no ID observations, year",g))
          }
          if(!is.finite(dpois(y.um[g,j,k],lamd.um[j]*K2D.sight[g,j,k],log=TRUE))){
            stop(paste("Starting observation model likelihood not finite. Unmarked observations, year",g))
          }
          if(!is.finite(dpois(y.unk[g,j,k],lamd.unk[j]*K2D.sight[g,j,k],log=TRUE))){
            stop(paste("Starting observation model likelihood not finite. Unknown marked status observations, year",g))
          }
        }
      }
      #integrated detector-session gamma-effect correction
      active.g <- z.super.init*z.init[,g]
      mu.det.g <- rep(0,J.sight[g])
      count.det.g <- rep(0,J.sight[g])
      for(k in 1:K.sight[g]){
        mu.det.g <- mu.det.g+
          K2D.sight[g,1:J.sight[g],k]*
          colSums(lamd[,1:J.sight[g],drop=FALSE]*active.g)
        count.det.g <- count.det.g+y.trap.total[g,1:J.sight[g],k]
      }
      logProb.trap.RE <- sum(lgamma(theta.d[g]+count.det.g)-lgamma(theta.d[g])+
                               theta.d[g]*log(theta.d[g])+mu.det.g-
                               (theta.d[g]+count.det.g)*log(theta.d[g]+mu.det.g))
      if(!is.finite(logProb.trap.RE)){
        stop(paste("Starting detector random-effect correction likelihood not finite, year",g))
      }
    }
  }
  dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are

  return(list(s=s.init,z=z.init,N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
              N.super=N.super.init,z.start=z.start.init,z.stop=z.stop.init,z.super=z.super.init,
              K1D.mark=K1D.mark,K2D.sight=K2D.sight,n.marked=n.marked,n.marked.all=n.marked.all,ID.marked=ID.marked,
              X.mark=X.mark,X.sight=X.sight,mark.states2D=mark.states2D,mark.states=mark.states,
              tel.z.states=tel.z.states,dummy.data=dummy.data,y2D=y2D,n.cap.all=n.cap.all,
              y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              y.trap.total=y.trap.total,trap.RE.dummy=trap.RE.dummy,
              xlim=xlim,ylim=ylim))
  
}
