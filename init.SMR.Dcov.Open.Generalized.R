e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.Open.Generalized <- function(data,inits=NA,M=NA){
  if(M < (data$n.marked.all)+1) stop("M must be larger than the number of marked individuals plus at least one unmarked individual.")
  library(abind)
  n.marked <- data$n.marked
  n.marked.all <- data$n.marked.all
  n.year <- data$n.year
  mark.states <- matrix(0,M,n.year)
  mark.states[1:n.marked.all,] <- data$mark.states
  #augment tel.z.states, code NA as 2
  tel.z.states <- matrix(2,M,n.year)
  tel.z.states[1:n.marked.all,] <- data$tel.z.states
  tel.z.states[is.na(tel.z.states)] <- 2  #use 2 to code NA for nimble function
  
  J.mark <- unlist(lapply(data$X.mark,nrow)) #traps per year
  J.sight <- unlist(lapply(data$X.sight,nrow)) #traps per year
  J.mark.max <- max(J.mark)
  J.sight.max <- max(J.sight)
  K.mark <- data$K.mark
  K.sight <- data$K.sight
  K.mark.max <- max(K.mark)
  K.sight.max <- max(K.sight)
  locs <- data$locs
  
  #augment y.mark and y.mnoID, pull out y.mnoID, y.um, y.unk
  y.mark <- array(0,dim=c(M,n.year,J.mark.max))
  y.mark[1:n.marked.all,,] <- data$y.mark
  y.mID <- array(0,dim=c(n.marked.all,n.year,J.sight.max))
  y.mID[1:n.marked.all,,] <- data$y.mID
  y.mnoID <- data$y.mnoID
  y.um <- data$y.um
  y.unk <- data$y.unk
  
  #reformat these
  ID.marked <- tel.inds <- matrix(0,max(n.marked),n.year)
  X.mark <- array(0,dim=c(n.year,J.mark.max,2))
  K1D.mark <- matrix(0,n.year,J.mark.max)
  X.sight <- array(0,dim=c(n.year,J.sight.max,2))
  K1D.sight <- matrix(0,n.year,J.sight.max)
  for(g in 1:n.year){
    if(n.marked[g]>0){
      ID.marked[1:n.marked[g],g] <- data$ID.marked[[g]]
      tel.inds[1:n.marked[g],g] <- data$tel.inds[[g]]
    }
    X.mark[g,1:J.mark[g],1:2] <- data$X.mark[[g]]
    K1D.mark[g,1:J.mark[g]] <- data$K1D.mark[[g]]
    X.sight[g,1:J.sight[g],1:2] <- data$X.sight[[g]]
    K1D.sight[g,1:J.sight[g]] <- data$K1D.sight[[g]]
  }
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  p0 <- inits$p0
  lam0 <- inits$lam0
  sigma <- inits$sigma
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  #but update s.inits for marked individuals before assigning latent detections
  idx <- which(rowSums(y.mID)>0)
  for(i in idx){
    trps <- matrix(0,nrow=0,ncol=2) #get locations of traps of capture across years for ind i
    for(g in 1:n.year){
      if(sum(y.mark[i,g,])>0){
        trps.g <- matrix(X.mark[g,which(y.mark[i,g,]>0),],ncol=2,byrow=FALSE)
        trps <- rbind(trps,trps.g)
      }
      if(sum(y.mID[i,g,])>0){
        trps.g <- matrix(X.sight[g,which(y.mID[i,g,]>0),],ncol=2,byrow=FALSE)
        trps <- rbind(trps,trps.g)
      }
    }
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  #update using telemetry if you have it
  if(!is.null(dim(data$locs))){
    n.tel.inds <- colSums(tel.inds>0)
    n.locs.ind <- matrix(0,max(n.tel.inds),n.year)
    for(g in 1:n.year){
      for(i in 1:n.tel.inds[g]){
        n.locs.ind[i,g]  <- sum(!is.na(locs[i,g,,1]))
      }
    }
    
    print("using telemetry to initialize telemetered s. Remove from data if not using in the model.")
    #update using telemetry if you have it
    tel.inds.all <- sort(unique(c(tel.inds)))
    idx <- which(tel.inds.all==0)
    if(length(idx)>0){
      tel.inds.all <- tel.inds.all[-idx]
    }
    
    for(i in tel.inds.all){
      locs.i <- matrix(0,nrow=0,ncol=2) #get locations of traps of capture across years for ind i
      for(g in 1:n.year){
        i.idx.this.g <- which(tel.inds[,g]==i)
        if(length(i.idx.this.g)>0){
          locs.g <- locs[i.idx.this.g,g,,]
          locs.i <- rbind(locs.i,locs.g)
        }
      }
      if(nrow(locs.i)>1){
        s.init[i,] <- colMeans(locs.i)
      }else{
        s.init[i,] <- locs.i
      }
      #make sure s is in state space
      if(s.init[i,1]<xlim[1]){
        s.init[i,1] <- xlim[1] + 0.01
      }
      if(s.init[i,1]>xlim[2]){
        s.init[i,1] <- xlim[2] - 0.01
      }
      if(s.init[i,2]<ylim[1]){
        s.init[i,2] <- ylim[1] + 0.01
      }
      if(s.init[i,2]>ylim[2]){
        s.init[i,2] <- ylim[2] - 0.01
      }
    }
  }else{
    tel.inds <- NA
    n.locs.ind <- NA
  }
  
  y.true <- array(0,dim=c(M,n.year,J.sight.max))
  y.true[1:n.marked.all,,] <- y.mID
  for(g in 1:n.year){
    D <- e2dist(s.init, X.sight[g,1:J.sight[g],])
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
      new.cell <- which(alldists[i,]==min(alldists[i,]))
      s.init[i,] <- data$dSS[new.cell,]
    }
  }
  
  #initialize z, start with possibly true detection history
  
  #initialize z, start with observed guys
  z.init <- matrix(0,M,n.year)
  z.start.init <- z.stop.init <- rep(0,M)
  y.true2D <- 1*(apply(y.true,c(1,2),sum)>0)
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
  z.obs <- 1*(rowSums(y.true)>0) #indicator for "ever observed"
  
  if(any(tel.z.states[z.init==1]==0))stop("At least one z initialized to 1 when tel.z.states=0. Bug in initialization code.")

  #initialize N structures from z.init
  N.init <- colSums(z.init[z.super.init==1,])
  N.survive.init <- N.recruit.init <- rep(NA,n.year-1)
  for(g in 2:n.year){
    N.survive.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1)
    N.recruit.init[g-1] <- N.init[g]-N.survive.init[g-1]
  }
  
  #get y2D constraints for z.start and z.stop update
  y.mark2D <- apply(y.mark,c(1,2),sum)
  y.mID2D <- apply(y.mID,c(1,2),sum)
  y2D <- apply(y.mID,c(1,2),sum) + apply(y.mark[1:n.marked.all,,],c(1,2),sum)
  #add telemetry states - using these instead of marked states since you can be marked and dead (how you observe telemetry death)
  for(i in 1:n.marked.all){
    y2D[data$tel.ID[i],which(data$tel.z.states[data$tel.ID[i],]>0)] <- 1
  }
  
  #check starting logProbs one session at a time
  for(g in 1:n.year){
    #marking process
    D.mark <- e2dist(s.init, X.mark[g,1:J.mark[g],])
    pd <- p0[g]*exp(-D.mark*D.mark/(2*sigma[g]*sigma[g]))
    logProb <- array(0,dim=c(M,J.mark[g]))
    for(i in 1:M){
      for(j in 1:J.mark[g]){
        logProb[i,j] <- dbinom(y.mark[i,g,j],size=K1D.mark[g,j],prob=pd[i,j],log=TRUE)
      }
    }
    if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marking process, year",g))
    #sighting process
    D.sight <- e2dist(s.init, X.sight[g,1:J.sight[g],])
    lamd <- lam0[g]*exp(-D.sight*D.sight/(2*sigma[g]*sigma[g]))
    #marked with ID obs
    if(n.marked[g]>0){
      logProb <- array(0,dim=c(n.marked[g],J.sight[g]))
      for(i in 1:n.marked[g]){
        for(j in 1:J.sight[g]){
          logProb[i,j] <- dpois(y.mID[i,g,j],lamd[i,j]*K1D.sight[g,j],log=TRUE)
        }
      }
      if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marked with ID observations, year",g))
    }
    #marked no ID obs
    logProb <- rep(0,J.sight[g])
    if(n.marked[g]>1){
      lamd.mnoID <- colSums(lamd[1:n.marked[g],1:J.sight[g]])
    }else if(n.marked[g]==1){
      lamd.mnoID <- lamd[n.marked[g],1:J.sight[g]]
    }else{#no marked guys
      lamd.mnoID <- rep(0,J.sight[g])
    }
    for(j in 1:J.sight[g]){
      logProb[j] <- dpois(y.mnoID[g,j],lamd.mnoID[j]*K1D.sight[g,j])
    }
    if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Marked with no ID observations, year",g))
    #um obs
    logProb <- rep(0,J.sight[g])
    lamd.um <- colSums(lamd[(n.marked[g]+1):M,1:J.sight[g]])
    for(j in 1:J.sight[g]){
      logProb[j] <- dpois(y.um[g,j],lamd.um[j]*K1D.sight[g,j])
    }
    if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Unmarked observations, year",g))
    #unk obs
    logProb <- rep(0,J.sight[g])
    lamd.unk <- colSums(lamd[1:M,1:J.sight[g]])
    for(j in 1:J.sight[g]){
      logProb[j] <- dpois(y.unk[g,j],lamd.unk[j]*K1D.sight[g,j])
    }
    if(!is.finite(sum(logProb)))stop(paste("Starting observation model likelihood not finite. Unknown marked status observations, year",g))
  }
  dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are

  
  return(list(s=s.init,z=z.init,N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
              N.super=N.super.init,z.start=z.start.init,z.stop=z.stop.init,z.super=z.super.init,
              K1D.mark=K1D.mark,K1D.sight=K1D.sight,n.marked=n.marked,n.marked.all=n.marked.all,ID.marked=ID.marked,
              tel.inds=tel.inds,X.mark=X.mark,X.sight=X.sight,mark.states=mark.states,tel.z.states=tel.z.states,
              dummy.data=dummy.data,y2D=y2D,
              y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              xlim=xlim,ylim=ylim,locs=locs,n.tel.inds=n.tel.inds,tel.inds=tel.inds,n.locs.ind=n.locs.ind))
  
}