e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.mobileAC.Open.Generalized.TrapRE <- function(data,inits=NA,M=NA){
  if(M < (data$n.cap.all)+1) stop("M must be larger than the number of captured individuals plus at least one unmarked individual.")
  library(abind)
  n.marked <- data$n.marked
  n.marked.all <- data$n.marked.all
  n.cap.all <- data$n.cap.all
  n.primary <- data$n.primary
  mark.states <- matrix(0,M,n.primary)
  mark.states[1:n.marked.all,] <- data$mark.states
  #augment tel.z.states, code NA as 2
  tel.z.states <- matrix(2,M,n.primary)
  tel.z.states[1:n.marked.all,] <- data$tel.z.states
  tel.z.states[is.na(tel.z.states)] <- 2  #use 2 to code NA for nimble function
  J.mark <- unlist(lapply(data$X.mark,nrow)) #traps per primary occasion
  J.sight <- unlist(lapply(data$X.sight,nrow)) #traps per primary occasion
  J.mark.max <- max(J.mark)
  J.sight.max <- max(J.sight)
  K.mark <- data$K.mark
  K.sight <- data$K.sight
  K.mark.max <- max(K.mark)
  K.sight.max <- max(K.sight)
  locs <- data$locs
  #augment y.mark and y.mnoID, pull out y.mnoID, y.um, y.unk
  y.mark <- array(0,dim=c(M,n.primary,J.mark.max))
  y.mark[1:n.cap.all,,] <- data$y.mark
  y.mID <- array(0,dim=c(n.marked.all,n.primary,J.sight.max))
  y.mID[1:n.marked.all,,] <- data$y.mID
  y.mnoID <- data$y.mnoID
  y.um <- data$y.um
  y.unk <- data$y.unk
  y.trap.total <- matrix(0,n.primary,J.sight.max)
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      y.trap.total[g,1:J.sight[g]] <- apply(y.mID[,g,1:J.sight[g],drop=FALSE],3,sum) +
        y.mnoID[g,1:J.sight[g]] +
        y.um[g,1:J.sight[g]] +
        y.unk[g,1:J.sight[g]]
    }
  }
  trap.RE.dummy <- rep(0,n.primary)
  
  #reformat these
  ID.marked <- matrix(0,max(n.marked),n.primary)
  X.mark <- array(0,dim=c(n.primary,J.mark.max,2))
  K1D.mark <- matrix(0,n.primary,J.mark.max)
  X.sight <- array(0,dim=c(n.primary,J.sight.max,2))
  K1D.sight <- matrix(0,n.primary,J.sight.max)
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
      K1D.sight[g,1:J.sight[g]] <- data$K1D.sight[[g]]
    }
  }
  xlim <- data$xlim
  ylim <- data$ylim
  ##pull out initial values
  p0 <- inits$p0
  lam0 <- inits$lam0
  sigma <- inits$sigma
  theta.d <- inits$theta.d
  sigma.move.init <- inits$sigma.move
  rsf.beta.init <- inits$rsf.beta
  D.beta1.init <- inits$D.beta1
  
  tau.move <- data$tau.move
  sigma.move.int <- sigma.move.init*sqrt(tau.move)
  
  e2dist <- function(x,y){
    i <- sort(rep(1:nrow(y),nrow(x)))
    dvec <- sqrt((x[,1]-y[i,1])^2 + (x[,2]-y[i,2])^2)
    matrix(dvec,nrow=nrow(x),ncol=nrow(y),byrow=FALSE)
  }
  n.cells <- nrow(data$dSS)
  n.cells.x <- length(data$x.vals)
  n.cells.y <- length(data$y.vals)
  #compute pi.cell from D.beta1.init
  lambda.cell <- data$InSS*exp(D.beta1.init*data$D.cov)
  pi.cell <- lambda.cell/sum(lambda.cell)
  #precompute rsf using rsf.beta.init
  rsf <- data$InSS*exp(rsf.beta.init*data$D.cov)
  y.true <- array(0,dim=c(M,n.primary,J.sight.max))
  y.true[1:n.marked.all,,] <- y.mID
  #provisional mobile activity centers for assigning latent detections
  z.super.pre <- rep(1,M)
  y.true.pre <- array(0,dim=c(M,n.primary,J.sight.max))
  y.true.pre[1:n.marked.all,,] <- y.mID
  s.pre <- array(0,dim=c(M,n.primary,2))
  avail.dist.pre <- use.dist.pre <- array(0,dim=c(M,n.primary-1,n.cells))
  on.inds <- which(z.super.pre==1)
  for(i in on.inds){
    obs2D <- rep(0,n.primary)
    for(g in 1:n.primary){
      obs2D[g] <- sum(y.mark[i,g,]) + sum(y.true.pre[i,g,])
    }
    if(!any(is.na(data$tel.ID))){
      tel.idx <- which(data$tel.ID==i)
      if(length(tel.idx)>0){
        for(tt in 1:data$n.tel.sessions[tel.idx]){
          gg <- data$tel.session[tel.idx,tt]
          if(!is.na(gg) && data$n.locs.ind[tel.idx,tt]>0){
            obs2D[gg] <- obs2D[gg] + data$n.locs.ind[tel.idx,tt]
          }
        }
      }
    }
    dets <- which(obs2D>0)
    if(length(dets)>0){
      first.det <- min(dets)
      last.det <- max(dets)
      #set detected primary occasions to mean trap/telemetry location
      for(g in dets){
        locs.g <- matrix(numeric(0),nrow=0,ncol=2)
        trapcaps <- which(y.mark[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.mark[[g]][trapcaps,,drop=FALSE])
        }
        trapcaps <- which(y.true.pre[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.sight[[g]][trapcaps,,drop=FALSE])
        }
        if(!any(is.na(data$tel.ID))){
          tel.idx <- which(data$tel.ID==i)
          if(length(tel.idx)>0){
            tel.g.idx <- which(data$tel.session[tel.idx,]==g)
            if(length(tel.g.idx)>0){
              nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
              if(nloc>0){
                locs.g <- rbind(locs.g,cbind(data$locs[tel.idx,tel.g.idx,1:nloc,1],data$locs[tel.idx,tel.g.idx,1:nloc,2]))
              }
            }
          }
        }
        mean.loc <- c(mean(locs.g[,1]),mean(locs.g[,2]))
        mean.loc[1] <- min(max(mean.loc[1],xlim[1]+1e-5),xlim[2]-1e-5)
        mean.loc[2] <- min(max(mean.loc[2],ylim[1]+1e-5),ylim[2]-1e-5)
        #check if mean location is in valid habitat
        mean.cell.x <- trunc(mean.loc[1]/data$res) + 1
        mean.cell.y <- trunc(mean.loc[2]/data$res) + 1
        mean.cell <- data$cells[mean.cell.x,mean.cell.y]
        if(data$InSS[mean.cell]==1){
          s.pre[i,g,] <- mean.loc
        }else{
          # snap to nearest InSS=1 cell by Euclidean distance
          dists <- sqrt((data$dSS[,1] - mean.loc[1])^2 + (data$dSS[,2] - mean.loc[2])^2)
          dists[data$InSS==0] <- Inf
          pick <- which.min(dists)
          s.pre[i,g,1] <- runif(1,data$dSS[pick,1] - data$res/2,data$dSS[pick,1] + data$res/2)
          s.pre[i,g,2] <- runif(1,data$dSS[pick,2] - data$res/2,data$dSS[pick,2] + data$res/2)
        }
      }
      #simulate backwards from first detection
      if(first.det>1){
        for(g in (first.det-1):1){
          avail <- getAvail(s=s.pre[i,g+1,],sigma=sigma.move.int[g],res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.pre[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.pre[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      #simulate forwards from last detection
      if(last.det<n.primary){
        for(g in (last.det+1):n.primary){
          avail <- getAvail(s=s.pre[i,g-1,],sigma=sigma.move.int[g-1],res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.pre[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.pre[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      # fill gaps with forward habitat-weighted simulation
      if(last.det>first.det){
        for(g in first.det:(last.det-1)){
          if(!(g+1)%in%dets){
            avail <- getAvail(s=s.pre[i,g,],sigma=sigma.move.int[g],res=data$res,
                              x.vals=data$x.vals,y.vals=data$y.vals,
                              n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
            use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
            s.cell <- sample(n.cells,1,prob=use)
            s.pre[i,g+1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
            s.pre[i,g+1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
          }
        }
      }
    }else{
      # undetected z.super=1: draw primary occasion 1 from pi.cell, subsequent from use distribution
      s.cell <- sample(n.cells,1,prob=pi.cell)
      s.pre[i,1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
      s.pre[i,1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      for(g in 2:n.primary){
        avail <- getAvail(s=s.pre[i,g-1,],sigma=sigma.move.int[g-1],res=data$res,
                          x.vals=data$x.vals,y.vals=data$y.vals,
                          n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
        use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
        s.cell <- sample(n.cells,1,prob=use)
        s.pre[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
        s.pre[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      }
    }
    for(g in 2:n.primary){
      avail.dist.pre[i,g-1,] <- getAvail(s=s.pre[i,g-1,1:2],sigma=sigma.move.int[g-1],res=data$res,
                                         x.vals=data$x.vals,y.vals=data$y.vals,
                                         n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
      use.dist.pre[i,g-1,] <- getUse(rsf=rsf,avail.dist=avail.dist.pre[i,g-1,],z.super=1)
    }
  }
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      D <- e2dist(s.pre[,g,],X.sight[g,1:J.sight[g],])
      lamd <- lam0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
      marked.inds <- which(mark.states[,g]==1)
      unmarked.inds <- which(mark.states[,g]==0)
      for(j in 1:J.sight[g]){
        #add marked no ID
        if(n.marked[g]>0){
          prob <- lamd[marked.inds,j]*(tel.z.states[marked.inds,g]!=0) #exclude known deaths
          prob <- prob/sum(prob)
          y.true[marked.inds,g,j] <- y.true[marked.inds,g,j] + rmultinom(1,y.mnoID[g,j],prob=prob)
        }
        #add unmarked
        prob <- lamd[,j]*(1-mark.states[,g])*(tel.z.states[,g]!=0) #zero out marked guys and dead guys
        prob <- prob/sum(prob)
        y.true[,g,j] <- y.true[,g,j] + rmultinom(1,y.um[g,j],prob=prob)
        #add unk
        prob <- lamd[,j]*(tel.z.states[,g]!=0) #exclude known deaths
        prob <- prob/sum(prob)
        y.true[,g,j] <- y.true[,g,j] + rmultinom(1,y.unk[g,j],prob=prob)
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
  #final mobile activity centers using allocated latent detections
  s.init <- array(0,dim=c(M,n.primary,2))
  avail.dist.init <- use.dist.init <- array(0,dim=c(M,n.primary-1,n.cells))
  on.inds <- which(z.super.init==1)
  for(i in on.inds){
    obs2D <- rep(0,n.primary)
    for(g in 1:n.primary){
      obs2D[g] <- sum(y.mark[i,g,]) + sum(y.true[i,g,])
    }
    if(!any(is.na(data$tel.ID))){
      tel.idx <- which(data$tel.ID==i)
      if(length(tel.idx)>0){
        for(tt in 1:data$n.tel.sessions[tel.idx]){
          gg <- data$tel.session[tel.idx,tt]
          if(!is.na(gg) && data$n.locs.ind[tel.idx,tt]>0){
            obs2D[gg] <- obs2D[gg] + data$n.locs.ind[tel.idx,tt]
          }
        }
      }
    }
    dets <- which(obs2D>0)
    if(length(dets)>0){
      first.det <- min(dets)
      last.det <- max(dets)
      #set detected primary occasions to mean trap/telemetry location
      for(g in dets){
        locs.g <- matrix(numeric(0),nrow=0,ncol=2)
        trapcaps <- which(y.mark[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.mark[[g]][trapcaps,,drop=FALSE])
        }
        trapcaps <- which(y.true[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.sight[[g]][trapcaps,,drop=FALSE])
        }
        if(!any(is.na(data$tel.ID))){
          tel.idx <- which(data$tel.ID==i)
          if(length(tel.idx)>0){
            tel.g.idx <- which(data$tel.session[tel.idx,]==g)
            if(length(tel.g.idx)>0){
              nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
              if(nloc>0){
                locs.g <- rbind(locs.g,cbind(data$locs[tel.idx,tel.g.idx,1:nloc,1],data$locs[tel.idx,tel.g.idx,1:nloc,2]))
              }
            }
          }
        }
        mean.loc <- c(mean(locs.g[,1]),mean(locs.g[,2]))
        mean.loc[1] <- min(max(mean.loc[1],xlim[1]+1e-5),xlim[2]-1e-5)
        mean.loc[2] <- min(max(mean.loc[2],ylim[1]+1e-5),ylim[2]-1e-5)
        #check if mean location is in valid habitat
        mean.cell.x <- trunc(mean.loc[1]/data$res) + 1
        mean.cell.y <- trunc(mean.loc[2]/data$res) + 1
        mean.cell <- data$cells[mean.cell.x,mean.cell.y]
        if(data$InSS[mean.cell]==1){
          s.init[i,g,] <- mean.loc
        }else{
          # snap to nearest InSS=1 cell by Euclidean distance
          dists <- sqrt((data$dSS[,1] - mean.loc[1])^2 + (data$dSS[,2] - mean.loc[2])^2)
          dists[data$InSS==0] <- Inf
          pick <- which.min(dists)
          s.init[i,g,1] <- runif(1,data$dSS[pick,1] - data$res/2,data$dSS[pick,1] + data$res/2)
          s.init[i,g,2] <- runif(1,data$dSS[pick,2] - data$res/2,data$dSS[pick,2] + data$res/2)
        }
      }
      #simulate backwards from first detection
      if(first.det>1){
        for(g in (first.det-1):1){
          avail <- getAvail(s=s.init[i,g+1,],sigma=sigma.move.int[g],res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.init[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.init[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      #simulate forwards from last detection
      if(last.det<n.primary){
        for(g in (last.det+1):n.primary){
          avail <- getAvail(s=s.init[i,g-1,],sigma=sigma.move.int[g-1],res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.init[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.init[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      # fill gaps with forward habitat-weighted simulation
      if(last.det>first.det){
        for(g in first.det:(last.det-1)){
          if(!(g+1)%in%dets){
            avail <- getAvail(s=s.init[i,g,],sigma=sigma.move.int[g],res=data$res,
                              x.vals=data$x.vals,y.vals=data$y.vals,
                              n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
            use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
            s.cell <- sample(n.cells,1,prob=use)
            s.init[i,g+1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
            s.init[i,g+1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
          }
        }
      }
    }else{
      # undetected z.super=1: draw primary occasion 1 from pi.cell, subsequent from use distribution
      s.cell <- sample(n.cells,1,prob=pi.cell)
      s.init[i,1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
      s.init[i,1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      for(g in 2:n.primary){
        avail <- getAvail(s=s.init[i,g-1,],sigma=sigma.move.int[g-1],res=data$res,
                          x.vals=data$x.vals,y.vals=data$y.vals,
                          n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
        use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
        s.cell <- sample(n.cells,1,prob=use)
        s.init[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
        s.init[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      }
    }
    for(g in 2:n.primary){
      avail.dist.init[i,g-1,] <- getAvail(s=s.init[i,g-1,1:2],sigma=sigma.move.int[g-1],res=data$res,
                                          x.vals=data$x.vals,y.vals=data$y.vals,
                                          n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
      use.dist.init[i,g-1,] <- getUse(rsf=rsf,avail.dist=avail.dist.init[i,g-1,],z.super=1)
    }
  }
  #check starting logProb
  logProb <- matrix(0,M,n.primary-1)
  for(i in 1:M){
    if(z.super.init[i]==1){
      for(g in 2:n.primary){
        logProb[i,g-1] <- dHabMove(x=s.init[i,g,1:2],s.prev=s.init[i,g-1,1:2],
                                   use.dist=use.dist.init[i,g-1,1:n.cells],
                                   dSS=data$dSS[1:n.cells,1:2],cells=data$cells[1:n.cells.x,1:n.cells.y],
                                   res=data$res,sigma.move=sigma.move.int[g-1],z.super=1,log=TRUE)
      }
    }
  }
  if(!all(is.finite(logProb))){ #can inspect LogProb object to identify problem individual-primary occasions
    stop("Starting logProb for activity centers is not finite, raise sigma.move.init. If that doesnt work, you may need to modify model or initialization algorithm.")
  }
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
  #add telemetry states - using these instead of marked states since you can be marked and dead (how you observe telemetry death)
  for(i in 1:n.marked.all){
    idx <- which(tel.z.states[i,]==1)
    if(length(idx)>0){
      y2D[i,idx] <- 1
    }
  }
  
  #z.start/z.stop sampler requires marked sighting primary occasions to fall
  #between the first and last primary occasions at which the individual is known alive
  for(i in 1:n.cap.all){
    dets <- which(y2D[i,]>0)
    if(length(dets)==0){
      stop(paste("Individual",i,"in 1:n.cap.all has no known-alive occasion in y2D."))
    }
    first.det <- min(dets)
    last.det <- max(dets)
    for(g in 1:n.primary){
      if(K.sight[g]>0){
        if((g<first.det | g>last.det) & mark.states[i,g]==1){
          stop(paste("Individual",i,"has mark.states==1 outside its known-alive interval;",
                     "this is incompatible with the z.start/z.stop sampler."))
        }
      }
    }
  }
  
  #check starting logProbs one session at a time
  for(g in 1:n.primary){
    #marking process
    if(J.mark[g]>0){
      D.mark <- e2dist(s.init[,g,],X.mark[g,1:J.mark[g],])
      pd <- p0[g]*exp(-D.mark*D.mark/(2*sigma[g]*sigma[g]))
      logProb <- array(0,dim=c(M,J.mark[g]))
      for(i in 1:M){
        for(j in 1:J.mark[g]){
          logProb[i,j] <- dbinom(y.mark[i,g,j],size=K1D.mark[g,j],prob=pd[i,j],log=TRUE)
        }
      }
      if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marking process, primary occasion",g))
    }
    #sighting process
    if(J.sight[g]>0){
      D.sight <- e2dist(s.init[,g,],X.sight[g,1:J.sight[g],])
      lamd <- lam0[g]*exp(-D.sight*D.sight/(2*sigma[g]*sigma[g]))
      #marked with ID obs
      if(n.marked[g]>0){
        marked.inds <- which(mark.states[,g]==1)
        logProb <- array(0,dim=c(length(marked.inds),J.sight[g]))
        for(ii in 1:length(marked.inds)){
          i <- marked.inds[ii]
          for(j in 1:J.sight[g]){
            logProb[ii,j] <- dpois(y.mID[i,g,j],lamd[i,j]*K1D.sight[g,j],log=TRUE)
          }
        }
        if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marked with ID observations, primary occasion",g))
      }
      #marked no ID obs
      logProb <- rep(0,J.sight[g])
      marked.inds <- which(mark.states[,g]==1)
      if(length(marked.inds)>1){
        lamd.mnoID <- colSums(lamd[marked.inds,1:J.sight[g],drop=FALSE])
      }else if(length(marked.inds)==1){
        lamd.mnoID <- lamd[marked.inds,1:J.sight[g]]
      }else{
        lamd.mnoID <- rep(0,J.sight[g])
      }
      for(j in 1:J.sight[g]){
        logProb[j] <- dpois(y.mnoID[g,j],lamd.mnoID[j]*K1D.sight[g,j])
      }
      if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marked with no ID observations, primary occasion",g))
      #um obs
      logProb <- rep(0,J.sight[g])
      unmarked.inds <- which(mark.states[,g]==0)
      lamd.um <- colSums(lamd[unmarked.inds,1:J.sight[g],drop=FALSE])
      for(j in 1:J.sight[g]){
        logProb[j] <- dpois(y.um[g,j],lamd.um[j]*K1D.sight[g,j])
      }
      if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Unmarked observations, primary occasion",g))
      #unk obs
      logProb <- rep(0,J.sight[g])
      lamd.unk <- colSums(lamd[1:M,1:J.sight[g]])
      for(j in 1:J.sight[g]){
        logProb[j] <- dpois(y.unk[g,j],lamd.unk[j]*K1D.sight[g,j])
      }
      if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Unknown marked status observations, primary occasion",g))
      #integrated detector-session gamma-effect correction
      active.g <- z.super.init*z.init[,g]
      mu.det.g <- K1D.sight[g,1:J.sight[g]]*colSums(lamd[,1:J.sight[g],drop=FALSE]*active.g)
      logProb.trap.RE <- sum(lgamma(theta.d[g] + y.trap.total[g,1:J.sight[g]]) - lgamma(theta.d[g]) + theta.d[g]*log(theta.d[g]) +
                               mu.det.g - (theta.d[g]+y.trap.total[g,1:J.sight[g]])*log(theta.d[g]+mu.det.g))
      if(!is.finite(logProb.trap.RE)){
        stop(paste("Starting detector random-effect correction likelihood not finite, primary occasion",g))
      }
    }
  }
  return(list(s=s.init,avail.dist=avail.dist.init,use.dist=use.dist.init,
              z=z.init,N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
              N.super=N.super.init,z.start=z.start.init,z.stop=z.stop.init,z.super=z.super.init,
              K1D.mark=K1D.mark,K1D.sight=K1D.sight,n.marked=n.marked,n.marked.all=n.marked.all,ID.marked=ID.marked,
              X.mark=X.mark,X.sight=X.sight,mark.states=mark.states,tel.z.states=tel.z.states,
              y2D=y2D,y.trap.total=y.trap.total,trap.RE.dummy=trap.RE.dummy,
              y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              xlim=xlim,ylim=ylim))
}