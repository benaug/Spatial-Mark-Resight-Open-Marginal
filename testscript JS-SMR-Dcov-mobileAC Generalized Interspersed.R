library(nimble)
library(coda)
source("sim.JS.SMR.Dcov.mobileAC.Generalized.Interspersed.R")
source("init.SMR.Dcov.mobileAC.Open.Generalized.Interspersed.R")
source("Nimble Model JS-SMR-Dcov-mobileAC Generalized Interspersed.R")
source("Nimble Functions JS-SMR-Dcov-mobileAC Generalized Interspersed.R") #contains custom distributions and updates
source("sSampler Dcov mobileAC Open Marginal Generalized Interspersed.R")
source("mask.check.R")
#must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

n.primary <- 6 #number of primary occasions
phi <- rep(0.8,n.primary-1) #per-capita recruitment by primary occasion
gamma <- rep(0.2,n.primary-1) #per-capita recruitment by primary occasion
tau <- rep(1,n.primary-1) #elapsed time between primary occasions for survival and recruitment
tau.move <- rep(1,n.primary-1) #movement time between primary occasions; set to 1 for discrete AC shifts
p0 <- rep(0.25,n.primary) #marking process p0
lam0 <- rep(0.25,n.primary) #sighting process lam0
sigma <- rep(0.5,n.primary) #detection function scale by primary occasion
sigma.move <- 2 #relocation scale by primary occasion
rsf.beta <- 0.5 #relocation RSF coefficient by primary occasion
p.mark <- rep(0.75,n.primary) #probability of marking given captured in marking process by primary occasion
#Number of occasions per primary occasion per method
#to skip sampling by a method in a primary occasion, set its K=0
K.mark <- c(8,8,8,8,8,8) #primary occasionly marking occasions
K.sight <- c(5,5,5,5,5,5) #primary occasionly resighting occasions
if(length(K.mark)!=length(K.sight))stop("K.mark and K.sight must be same length")
if(length(K.mark)!=n.primary)stop("K.mark and K.sight must be of length n.primary")

#theta is probability of observing each sample type for marked and unmarked individuals
#assuming the same over primary occasions
theta.marked <- c(0.95,0.025,0.025) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.95 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)

#make an SCR trapping array. Making the trapping array size vary by session
#For occasions with no marking or sighting, insert a trap matrix with 0 rows, matrix(0,nrow=0,ncol=2)
X.sight <- vector("list",n.primary)
X.sight[[1]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[2]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[3]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[4]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[5]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[6]] <- as.matrix(expand.grid(1:10,1:10))

#Let's spread the marking traps out across the study area by randomly selecting from the 
# sighting traps. For movement, this is probably better than only allowing marked and thus identifiable
#individuals in a central subregion of the state space only, like other data simulators in this repo are set up for
X.mark <- vector("list",n.primary)
X.mark[[1]] <- X.sight[[1]][sort(sample(1:nrow(X.sight[[1]]),36)),]
X.mark[[2]] <- X.sight[[2]][sort(sample(1:nrow(X.sight[[2]]),36)),]
X.mark[[3]] <- X.sight[[3]][sort(sample(1:nrow(X.sight[[3]]),36)),]
X.mark[[4]] <- X.sight[[4]][sort(sample(1:nrow(X.sight[[4]]),36)),]
X.mark[[5]] <- X.sight[[5]][sort(sample(1:nrow(X.sight[[5]]),36)),]
X.mark[[6]] <- X.sight[[6]][sort(sample(1:nrow(X.sight[[6]]),36)),]

#Check for consistency between traps and occasions
for(g in 1:n.primary){
  if(K.sight[g]==0){ #trap matrix should have 0 rows
    if(!nrow(X.sight[[g]])==0){
      stop(paste("X.sight and K.sight inconsistent, session",g))
    }
  }else{ #trap matrix should have >0 rows
    if(nrow(X.sight[[g]])==0){
      stop(paste("X.sight and K.sight inconsistent, session",g))
    }
  }
  if(K.mark[g]==0){ #trap matrix should have 0 rows
    if(!nrow(X.mark[[g]])==0){
      stop(paste("X.mark and K.mark inconsistent, session",g))
    }
  }else{ #trap matrix should have >0 rows
    if(nrow(X.mark[[g]])==0){
      stop(paste("X.mark and K.mark inconsistent, session",g))
    }
  }
}

#For each primary session, a character vector of length K.mark[g] + K.sight[g] specifying the order of the marking
#and sighting occasions. Vector elements are either "M" or "S" arranged in the order the marking
#and sighting occasions occurred. There must be K.mark[g] M's and K.sight[g] S's. 
#Set K.order[[g]] <- NA only if there is no marking or sighting effort in primary occasion g
#Data simulator will check these requirements.
K.order <- vector("list",n.primary)
K.order[[1]] <- c("M","M","M","M","M","S","S","S","M","M","M","S","S")
K.order[[2]] <-  c("M","M","M","M","M","S","S","S","M","M","M","S","S")
K.order[[3]] <-  c("M","M","M","M","M","S","S","S","M","M","M","S","S")
K.order[[4]] <- c("M","M","M","M","S","S","S","M","M","M","M","S","S")
K.order[[5]] <- c("M","M","M","M","S","S","S","M","M","M","M","S","S")
K.order[[6]] <-c("M","M","M","M","S","S","S","M","M","M","M","S","S")


### Habitat covariate stuff###
#get x and y extent for each grid separately, then merge
xlim <- ylim <- matrix(NA,n.primary,2)
buff <- 3 #state space buffer around traps
X.both <- vector("list",n.primary)
for(g in 1:n.primary){
  X.both[[g]] <- rbind(X.mark[[g]],X.sight[[g]])
  if(nrow(X.both[[g]])>0){
    xlim[g,] <- range(X.both[[g]][,1]) + c(-buff,buff)
    ylim[g,] <- range(X.both[[g]][,2]) + c(-buff,buff)
  }
}
xlim <- c(min(xlim[,1],na.rm=TRUE),max(xlim[,2],na.rm=TRUE))
ylim <- c(min(ylim[,1],na.rm=TRUE),max(ylim[,2],na.rm=TRUE))

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim - x.shift
ylim <- ylim - y.shift
for(g in 1:n.primary){
  X.both[[g]][,1] <- X.both[[g]][,1] - x.shift
  X.both[[g]][,2] <- X.both[[g]][,2] - y.shift
  X.mark[[g]][,1] <- X.mark[[g]][,1] - x.shift
  X.mark[[g]][,2] <- X.mark[[g]][,2] - y.shift
  X.sight[[g]][,1] <- X.sight[[g]][,1] - x.shift
  X.sight[[g]][,2] <- X.sight[[g]][,2] - y.shift
}

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#for plotting, making mask
X.mark.all <- X.sight.all <- matrix(NA,nrow=0,ncol=2)
for(g in 1:n.primary){
  X.mark.all <- rbind(X.mark.all,X.mark[[g]])
  X.sight.all <- rbind(X.sight.all,X.sight[[g]])
}
X.all <- rbind(X.mark.all,X.sight.all)

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
set.seed(1333)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(5,5),messages=FALSE)[[2]]
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X.sight.all,pch=4,lwd=2)
points(X.mark.all,pch=4,lwd=2,col="darkred")

#Additionally, maybe we want to exclude "non-habitat" or limit the state space extent
#let's use a 3sigma buffer
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(0,length(D.cov))
dists <- e2dist(X.all,dSS.tmp)
min.dists <- apply(dists,2,min)
InSS[min.dists<(3*max(sigma))] <- 1
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="Habitat",col=cols1)
points(X.all,pch=4,col="darkred",lwd=2)

#Density covariates
D.beta0 <- -0.75 #data simulator uses intercept for marked + unmarked
D.beta1 <- 0.5
#what is implied expected primary occasion 1 N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected primary occasion 1 N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density in primary occasion 1",col=cols1)
points(X.sight.all,pch=4,cex=1,lwd=2)
points(X.mark.all,pch=4,cex=1,lwd=2,col="darkred")

#Mark/Telemetry settings - For now, we assume mark history is known and deaths observed if marked at time of death
#this is a simplified scenario for using telemetry collars for marks
n.tel.locs <- 15 #number of locs per individual
mark.g.pars <- c(2,2,3) #parameters for truncated poisson: c(lambda, lower truncation, upper truncation)
#data simulator requires lower bound be 1 or higher. 1 means it fails before 2nd primary occasion
#mark lifetime frequencies for mark.g.pars
table(rtruncpois(10000,lambda=mark.g.pars[1],lower=mark.g.pars[2],upper=mark.g.pars[3]))/10000
#marking protocol: #1) never replace a mark if currently collared on capture 2) always replace
mark.protocol <- 2 

# simulate some data
set.seed(390297) #change seed for new data set
data <- sim.JS.SMR.Dcov.mobileAC.Generalized.Interspersed(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,
                                                          tau=tau,tau.move=tau.move,
                                                          InSS=InSS,phi=phi,gamma=gamma,n.primary=n.primary,K.order=K.order,
                                                          theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                                                          p0=p0,lam0=lam0,sigma=sigma,sigma.move=sigma.move,rsf.beta=rsf.beta,
                                                          K.mark=K.mark,K.sight=K.sight,
                                                          X.mark=X.mark,X.sight=X.sight,xlim=xlim,ylim=ylim,res=res,
                                                          mark.g.pars=mark.g.pars,mark.protocol=mark.protocol,
                                                          p.mark=p.mark,n.tel.locs=n.tel.locs)

#what is observed data? Note data objects have all n.primarys with all 0 data if no effort for a method
#Could be structured without primary occasions with no effort, but that would require more work changing custom
#N/z updates.

#mark and sight data summed over occasions
#str(data$y.mark) #marking process history: n.cap.all x n.primary x J.mark.max.
#total number captured (n.cap.all) might be > total number ever marked (n.marked.all). 
#if so, marked individuals must be first, then captured but unmarked individuals
#str(data$y.mID) #marked with ID sighting history: n.marked.all x n.primary x J.sight.max
#str(data$y.mnoID) #marked with no ID sighting history: n.primary x J.sight.max
#str(data$y.um) #unmarked sighting history: n.primary x J.sight.max
#str(data$y.unk) #unknown marked status sighting history: n.primary x J.sight.max
#str(data$mark.states) #mark status history: n.marked.all x n.primary
#str(data$tel.z.states) #telemetry survival observations: n.marked.all x n.primary
#str(data$locs) #telemetry locations: n.tel.inds x n.tel.sessions.max x n.tel.locs.max x 2
#use tel.ID and tel.session to map to individual and population primary occasion


#these plots are cool, you should look at these.

#visualize expected relative density and realized activity centers in each primary occasion
#primary occasion 1: cell colors depict expected relative density in primary occasion 1
#primary occasions 2 on: cell colors depict expected relative density given the realized s[g-1] and z[g-1] from previous primary occasion
#points are realized activity centers for this expectation
# par(mfrow=c(1,1),ask=FALSE)
# for(plot.year in 1:n.primary){
#   image(x.vals,y.vals,matrix(data$truth$pi.cell[plot.year,],n.cells.x,n.cells.y),
#         main=paste("Expected Relative Density, Year", plot.year),
#         col=cols1)
#   points(X.all,pch=4,cex=0.75)
#   points(data$truth$s[data$truth$z[,plot.year]==1,plot.year,],pch=16)
# }

#visualize individual movement trajectories. Start primary occasion is larger circle
#will be a mess with a lot of individuals.
# ind.cols <- c("#E63946","#FF9F1C","#FFDD00","#2EC4B6","#3A86FF",
#               "#8338EC","#FB5607","#06D6A0","#FFB700","#118AB2",
#               "#EF476F","#FFC8DD","#B5E48C","#00BBF9","#9B5DE5",
#               "#F15BB5","#00F5D4","#FEE440","#FF595E","#6A4C93")
# n.colors <- length(ind.cols)
# par(mfrow=c(1,1), ask=FALSE)
# #can use just InSS to make sure all s inside state space
# image(x.vals, y.vals, matrix(D.cov*InSS, n.cells.x, n.cells.y), col=cols1,
#       main="Movement Trajectories over D.cov")
# for(i in 1:data$truth$N.super){
#   #skip individuals alive for less than 2 primary occasions
#   if(sum(data$truth$z[i,])<2) next
#   ind.col <- ind.cols[(i-1) %% n.colors + 1]
#   alive.years <- which(data$truth$z[i,]==1) #get primary occasions alive
#   # plot points for all alive primary occasions
#   points(data$truth$s[i,alive.years,1],
#          data$truth$s[i,alive.years,2],
#          pch=16, col=ind.col, cex=0.8)
# 
#   #plot lines between consecutive alive primary occasions
#   if(length(alive.years)>1){
#     for(t in 1:(length(alive.years)-1)){
#       if(alive.years[t+1] == alive.years[t]+1){
#         lines(x=c(data$truth$s[i,alive.years[t],1],data$truth$s[i,alive.years[t+1],1]),
#               y=c(data$truth$s[i,alive.years[t],2],data$truth$s[i,alive.years[t+1],2]),
#               col=ind.col,lwd=1.2)
#       }
#     }
#   }
#   #mark entry point with larger circle
#   points(data$truth$s[i,alive.years[1],1],data$truth$s[i,alive.years[1],2],
#          pch=16,col=ind.col,cex=1.5)
# }
# points(X.all,pch=4,cex=0.75,lwd=2)

# can look at individual by primary occasion availability and use distributions
# i <- 1
# g <- 1
# par(mfrow=c(3,1))
# image(x.vals,y.vals,matrix(data$truth$avail.dist[i,g,],n.cells.x,n.cells.y),
#       main="Availability Distribution",col=cols1)
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="darkred")
# image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="RSF Cov",col=cols1)
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="darkred")
# image(x.vals,y.vals,matrix(data$truth$use.dist[i,g,],n.cells.x,n.cells.y),col=cols1,
#       main="Use Distribution")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="darkred")
# par(mfrow=c(1,1))


data$N #abundance by primary occasion
colSums(apply(data$y.mark>0,c(1,2),sum)>0) #total marking process captures per primary occasion
colSums(data$mark.deploy) #total marks deployed per primary occasion
rowSums(data$mark.deploy) #total marks deployed per captured individual
data$n.marked #marks active per primary occasion

#total detected individuals
colSums(apply(data$truth$y,c(1,2),sum)>0)
#marked spatial recaps
table(apply(1*(data$truth$y.mark>0),c(1,2),sum))

#visualize detections by primary occasion. only showing SCR and identified SMR detections
for(g in 1:n.primary){
  image(data$x.vals,data$y.vals,matrix(data$D.cov*data$InSS,data$n.cells.x,data$n.cells.y),
        main=paste("primary occasion",g),xlab="X",ylab="Y",col=cols1)
  if(data$J.sight[g]>0){
    points(data$X.sight[[g]],pch=4,lwd=2)
  }
  if(data$J.mark[g]>0){
    points(data$X.mark[[g]],pch=4,lwd=2,col="darkred")
  }
  alive.g <- which(data$truth$z[,g]==1)
  points(data$truth$s[alive.g,g,1],data$truth$s[alive.g,g,2],pch=16)
  if(data$n.marked[g]>0){
    for(i in 1:data$n.marked[g]){
      id <- data$ID.marked[[g]][i]
      traps <- matrix(numeric(0),nrow=0,ncol=2)
      if(data$J.sight[g]>0){
        trapcaps <- which(rowSums(data$y.mID[id,g,,])>0)
        if(length(trapcaps)>0){
          traps <- rbind(traps,data$X.sight[[g]][1:data$J.sight[g],][trapcaps,])
        }
      }
      if(data$J.mark[g]>0){
        trapcaps2 <- which(data$y.mark[id,g,]>0)
        if(length(trapcaps2)>0){
          traps <- rbind(traps,data$X.mark[[g]][1:data$J.mark[g],][trapcaps2,])
        }
      }
      s <- data$s[id,g,]
      points(s[1],s[2],col="goldenrod",pch=16)
      if(nrow(traps)>0){
        for(j in 1:nrow(traps)){
          lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
        }
      }
      tel.idx <- which(data$tel.ID==id)
      if(length(tel.idx)>0){
        tel.g.idx <- which(data$tel.session[tel.idx,]==g)
        if(length(tel.g.idx)>0){
          nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
          if(nloc>0){
            for(l in 1:nloc){
              lines(x=c(s[1],data$locs[tel.idx,tel.g.idx,l,1]),
                    y=c(s[2],data$locs[tel.idx,tel.g.idx,l,2]),
                    col="gray80")
            }
            points(data$locs[tel.idx,tel.g.idx,1:nloc,1],data$locs[tel.idx,tel.g.idx,1:nloc,2],
                   pch=16,cex=0.5,col="lightblue")
            points(s[1],s[2],col="darkblue",pch=16)
          }
        }
      }
    }
  }
}

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

##Initialize##
data$N[1] + sum(data$N.recruit) #true N.super

M <- 250 #data augmentation level.

#initialize N and z objects and activity centers
if(M < (data$n.marked.all)+1) stop("M must be larger than the number of marked individuals plus at least one unmarked individual.")
#pull these from data (won't be in environment if not simulated directly above)
n.primary <- data$n.primary #number of primary sessions
n.marked <- data$n.marked #number of individuals carrying a mark in each primary occasion
n.marked.all <- data$n.marked.all #total number of individuals ever marked
n.cap.all <- data$n.cap.all #total number of individuals ever captured (might be every marked individual)
J.mark <- data$J.mark
J.sight <- data$J.sight
K.mark <- data$K.mark
K.sight <- data$K.sight
xlim <- data$xlim
ylim <- data$ylim
dSS <- data$dSS
cells <- data$cells
res <- data$res
cellArea <- res^2
D.cov <- data$D.cov
InSS <- data$InSS
x.vals <- data$x.vals
y.vals <- data$y.vals
n.cells <- data$n.cells
n.cells.x <- data$n.cells.x
n.cells.y <- data$n.cells.y
max.n.tel.locs <- dim(data$locs)[3]
#make sure any NA's in locs are converted to 0 to avoid nimble warnings about NA data
idx <- which(is.na(data$locs))
if(length(idx)>0){
  data$locs[idx] <- 0
}
tau <- data$tau
tau.move <- data$tau.move

#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
inits <- list(p0=rep(0.1,n.primary),lam0=rep(0.25,n.primary),#initializing with 1 parameter per session, just set all to same value
              sigma=rep(0.5,n.primary),
              sigma.move=2,D.beta1=0.5,rsf.beta=0.5) #single sigma.move, D.beta1, rsf.beta
#This function structures the simulated data to fit the model in Nimble (some more restructing below)
nimbuild <- init.SMR.Dcov.mobileAC.Open.Generalized.Interspersed(data,inits,M=M)

#plot to check s inits by primary occasion
for(g in 1:n.primary){
  image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),
        main=paste("primary occasion",g),xlab="X",ylab="Y",col=cols1)
  if(J.sight[g]>0){
    points(data$X.sight[[g]],pch=4,lwd=2)
  }
  if(J.mark[g]>0){
    points(data$X.mark[[g]],pch=4,lwd=2,col="darkred")
  }
  alive.g <- which(nimbuild$z[,g]==1)
  points(nimbuild$s[alive.g,g,1],nimbuild$s[alive.g,g,2],pch=16) #initialized activity centers
  if(n.marked[g]>0){
    for(i in 1:n.marked[g]){
      id <- data$ID.marked[[g]][i]
      traps <- matrix(numeric(0),nrow=0,ncol=2)
      if(data$J.sight[g]>0){
        trapcaps <- which(rowSums(data$y.mID[id,g,,])>0)
        if(length(trapcaps)>0){
          traps <- rbind(traps,data$X.sight[[g]][1:data$J.sight[g],][trapcaps,])
        }
      }
      if(data$J.mark[g]>0){
        trapcaps2 <- which(data$y.mark[id,g,]>0)
        if(length(trapcaps2)>0){
          traps <- rbind(traps,data$X.mark[[g]][1:data$J.mark[g],][trapcaps2,])
        }
      }
      s <- nimbuild$s[id,g,1:2]
      points(s[1],s[2],col="goldenrod",pch=16)
      if(nrow(traps)>0){
        for(j in 1:nrow(traps)){
          lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
        }
      }
    }
  }
}

#these indicate in which year marking/sighting occurs and how many total sessions of each
mark.g <- which(K.mark!=0)
sight.g <- which(K.sight!=0)
n.mark.g <- length(mark.g)
n.sight.g <- length(sight.g)

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.primary=n.primary,M=M,J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,
                  K1D.mark=nimbuild$K1D.mark,K2D.sight=nimbuild$K2D.sight,
                  D.cov=D.cov,tau=tau,tau.move=tau.move,
                  n.marked.all=nimbuild$n.marked.all,
                  n.tel.sessions=data$n.tel.sessions,tel.session=data$tel.session,max.n.tel.locs=max.n.tel.locs,
                  tel.ID=data$tel.ID,n.tel.inds=data$n.tel.inds,n.locs.ind=data$n.locs.ind,
                  mark.g=mark.g,sight.g=sight.g,n.mark.g=n.mark.g,
                  n.sight.g=n.sight.g,n.cells=n.cells,n.cells.x=n.cells.x,
                  n.cells.y=n.cells.y,res=res,x.vals=x.vals,y.vals=y.vals,xlim=xlim,ylim=ylim,
                  cellArea=cellArea)
#inits for Nimble
Niminits <- list(N=nimbuild$N,N.survive=nimbuild$N.survive,N.recruit=nimbuild$N.recruit,
                 ER=nimbuild$N.recruit,N.super=nimbuild$N.super,z.super=nimbuild$z.super,
                 z=nimbuild$z,z.start=nimbuild$z.start,z.stop=nimbuild$z.stop,
                 s=nimbuild$s,phi.fixed=0.5,D0=nimbuild$N[1]/(sum(InSS)*res^2),
                 p0=inits$p0[mark.g],lam0=inits$lam0[sight.g],sigma.fixed=inits$sigma[1],
                 sigma.move=inits$sigma.move,rsf.beta=inits$rsf.beta,D.beta1=inits$D.beta1)

#data for Nimble
Nimdata <- list(y.mark=nimbuild$y.mark, #marking process
                y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                mark.states=nimbuild$mark.states, #mark state history (who is marked in each year)
                tel.z.states=nimbuild$tel.z.states, #telemetry z state observations
                cells=cells,InSS=InSS,dSS=dSS,
                X.mark=nimbuild$X.mark,X.sight=nimbuild$X.sight,locs=data$locs)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','N.super','lambda.y1',
                'phi.fixed','p0','lam0','sigma.fixed','theta.marked','theta.unmarked',
                'D0','D.beta1','sigma.move','rsf.beta')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('phi.fixed','gamma','p0','lam0','sigma.fixed',
                  'theta.marked','theta.unmarked[2:3]','sigma.move','rsf.beta')
#use this above if theta.unmarked is primary occasion specific (no change for theta.marked): paste('theta.unmarked[1:',n.primary,',2:3',']')
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,nodes=config.nodes)

#add N/z sampler
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration?
#25% of M seems reasonable, but optimal will depend on data set

y.mark.nodes <- Rmodel$expandNodeNames(paste0("y.mark[1:",M,",1:",n.primary,",1:",max(J.mark),"]"))
pd.nodes <- Rmodel$expandNodeNames(paste0("pd[1:",M,",1:",n.primary,",1:",max(J.mark),"]"))
lam.nodes <- Rmodel$expandNodeNames(paste0("lam[1:",M,",1:",n.primary,",1:",max(J.sight),"]"))

# One vector node for each sighting occasion k; build explicitly so node order is
# known for the partial z.start/z.stop likelihood recalculations.
lam.um.nodes <- lam.unk.nodes <- y.um.nodes <- y.unk.nodes <- c()
sight.node.start <- integer(n.sight.g)
for(g2 in 1:n.sight.g){
  gg <- sight.g[g2]
  sight.node.start[g2] <- length(y.um.nodes)+1L
  for(k in 1:K.sight[gg]){
    y.um.tmp <- Rmodel$expandNodeNames(paste0("y.um[",gg,",1:",J.sight[gg],",",k,"]"))
    y.unk.tmp <- Rmodel$expandNodeNames(paste0("y.unk[",gg,",1:",J.sight[gg],",",k,"]"))
    lam.um.tmp <- Rmodel$expandNodeNames(paste0("lam.um[",gg,",1:",J.sight[gg],",",k,"]"))
    lam.unk.tmp <- Rmodel$expandNodeNames(paste0("lam.unk[",gg,",1:",J.sight[gg],",",k,"]"))
    if(length(y.um.tmp)!=1|length(y.unk.tmp)!=1|length(lam.um.tmp)!=1|length(lam.unk.tmp)!=1){
      stop("Expected one vector node per sighting occasion for y.um, y.unk, lam.um, and lam.unk.")
    }
    y.um.nodes <- c(y.um.nodes,y.um.tmp)
    y.unk.nodes <- c(y.unk.nodes,y.unk.tmp)
    lam.um.nodes <- c(lam.um.nodes,lam.um.tmp)
    lam.unk.nodes <- c(lam.unk.nodes,lam.unk.tmp)
  }
}

N.nodes <- Rmodel$expandNodeNames("N")
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.primary-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.primary-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.primary-1,"]"))
s.nodes <- Rmodel$expandNodeNames("s")
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
tel.z.states.nodes <- Rmodel$expandNodeNames(paste0("tel.z.states[1:",M,",1]"))
calcNodes <- c(s.nodes,N.nodes,N.recruit.nodes,y.mark.nodes,y.um.nodes,y.unk.nodes,z.nodes,tel.z.states.nodes)

# cells must be double inside compiled custom sampler
cells.double <- matrix(as.double(cells),n.cells.x,n.cells.y)

conf$addSampler(target=c("z"),
                type='zSampler',
                control=list(M=M,n.marked.all=n.marked.all,n.cap.all=n.cap.all,
                             n.primary=n.primary,J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,
                             mark.g=mark.g,sight.g=sight.g,n.mark.g=n.mark.g,n.sight.g=n.sight.g,
                             xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,
                             cells=cells.double,dSS=dSS,n.cells=n.cells,
                             n.cells.x=n.cells.x,n.cells.y=n.cells.y,res=res,
                             mark.states=nimbuild$mark.states,
                             tel.z.states=nimbuild$tel.z.states,
                             z.super.ups=z.super.ups,y2D=nimbuild$y2D,
                             y.mark.nodes=y.mark.nodes,pd.nodes=pd.nodes,
                             y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                             lam.nodes=lam.nodes,tel.z.states.nodes=tel.z.states.nodes,
                             lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                             sight.node.start=sight.node.start,
                             N.nodes=N.nodes,z.nodes=z.nodes,ER.nodes=ER.nodes,
                             N.survive.nodes=N.survive.nodes,N.recruit.nodes=N.recruit.nodes,
                             s.nodes=s.nodes,calcNodes=calcNodes),silent=TRUE)


#activity center samplers. There are 2 samplers here for these cases
#1) z.super=1 and z=1, sSampler1 uses Metropolis-Hastings
#2) z.super=1 and z=0, sSampler2 uses Metropolis-Hastings with proposal sd tuned separately from above
#z.super=0, do nothing, activity centers likelihood is gated and all values are set to 0
# Both implement Herliansyah et al. (2024) efficiency improvement
#that does not require summation over all individuals' lam values when updating 1 activity center.
for(i in 1:M){
  for(g in 1:n.primary){
    s.target <- paste0("s[",i,",",g,",1:2]")
    calcNodes <- Rmodel$getDependencies(s.target)
    s.nodes <- Rmodel$expandNodeNames(s.target)
    if(g<n.primary){
      s.nodes <- c(s.nodes,
                   Rmodel$expandNodeNames(paste0("avail.dist[",i,",",g,",1:",n.cells,"]")),
                   Rmodel$expandNodeNames(paste0("s[",i,",",g+1,",1:2]")))
    }
    loc.nodes <- c()
    tel.idx <- which(data$tel.ID==i)
    if(length(tel.idx)>0){
      tel.g.idx <- which(data$tel.session[tel.idx,]==g)
      if(length(tel.g.idx)>0){
        nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
        if(nloc>0){
          loc.nodes <- c(loc.nodes,Rmodel$expandNodeNames(paste0("locs[",tel.idx,",",tel.g.idx,",1:",nloc,",1:2]")))
        }
      }
    }
    if(length(loc.nodes)>0){
      s.nodes <- c(s.nodes,loc.nodes)
    }
    s.nodes <- unique(s.nodes)
    if(J.mark[g]>0){
      pd.nodes <- Rmodel$expandNodeNames(paste0("pd[",i,",",g,",1:",J.mark[g],"]"))
      y.mark.nodes <- Rmodel$expandNodeNames(paste0("y.mark[",i,",",g,",1:",J.mark[g],"]"))
    }else{
      pd.nodes <- character(0)
      y.mark.nodes <- character(0)
    }
    if(J.sight[g]>0 & K.sight[g]>0){
      lam.nodes <- Rmodel$expandNodeNames(paste0("lam[",i,",",g,",1:",J.sight[g],"]"))
      y.um.nodes <- Rmodel$expandNodeNames(paste0("y.um[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
      y.unk.nodes <- Rmodel$expandNodeNames(paste0("y.unk[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
      lam.um.nodes <- Rmodel$expandNodeNames(paste0("lam.um[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
      lam.unk.nodes <- Rmodel$expandNodeNames(paste0("lam.unk[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
      if(i<=nimbuild$n.marked.all){
        y.mID.nodes <- Rmodel$expandNodeNames(paste0("y.mID[",i,",",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
        y.mnoID.nodes <- Rmodel$expandNodeNames(paste0("y.mnoID[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
        lam.mnoID.nodes <- Rmodel$expandNodeNames(paste0("lam.mnoID[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
      }else{
        y.mID.nodes <- character(0)
        y.mnoID.nodes <- character(0)
        lam.mnoID.nodes <- character(0)
      }
    }else{
      lam.nodes <- character(0)
      y.mID.nodes <- character(0)
      y.mnoID.nodes <- character(0)
      y.um.nodes <- character(0)
      y.unk.nodes <- character(0)
      lam.mnoID.nodes <- character(0)
      lam.um.nodes <- character(0)
      lam.unk.nodes <- character(0)
    }
    conf$addSampler(target=s.target,
                    type='sSampler1',
                    control=list(i=i,g=g,xlim=xlim,ylim=ylim,J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,
                                 n.marked.all=nimbuild$n.marked.all,mark.states=nimbuild$mark.states[i,,],
                                 s.nodes=s.nodes,pd.nodes=pd.nodes,mark.states.all=nimbuild$mark.states,
                                 lam.nodes=lam.nodes,y.mark.nodes=y.mark.nodes,
                                 y.mID.nodes=y.mID.nodes,y.mnoID.nodes=y.mnoID.nodes,
                                 y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                 lam.mnoID.nodes=lam.mnoID.nodes,lam.um.nodes=lam.um.nodes,
                                 lam.unk.nodes=lam.unk.nodes,calcNodes=calcNodes,scale=1),
                    silent=TRUE)
    conf$addSampler(target=s.target,
                    type='sSampler2',
                    control=list(i=i,g=g,xlim=xlim,ylim=ylim,
                                 s.nodes=s.nodes,scale=1),
                    silent=TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#optional gap-jumper activity center sampler
#3) z.super=1, either z state, "gap jumper" to better jump across state space gaps
#proposal sigma is sigma.move*jump.multiplier so it scales with sigma.move.
# for(i in 1:M){
#   for(g in 1:n.primary){
#     s.target <- paste0("s[",i,",",g,",1:2]")
#     calcNodes <- Rmodel$getDependencies(s.target)
#     s.nodes <- Rmodel$expandNodeNames(s.target)
#     if(g<n.primary){
#       s.nodes <- c(s.nodes,
#                    Rmodel$expandNodeNames(paste0("avail.dist[",i,",",g,",1:",n.cells,"]")),
#                    Rmodel$expandNodeNames(paste0("s[",i,",",g+1,",1:2]")))
#     }
#     loc.nodes <- c()
#     tel.idx <- which(data$tel.ID==i)
#     if(length(tel.idx)>0){
#       tel.g.idx <- which(data$tel.session[tel.idx,]==g)
#       if(length(tel.g.idx)>0){
#         nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
#         if(nloc>0){
#           loc.nodes <- c(loc.nodes,Rmodel$expandNodeNames(paste0("locs[",tel.idx,",",tel.g.idx,",1:",nloc,",1:2]")))
#         }
#       }
#     }
#     if(length(loc.nodes)>0){
#       s.nodes <- c(s.nodes,loc.nodes)
#     }
#     s.nodes <- unique(s.nodes)
#     if(J.mark[g]>0){
#       pd.nodes <- Rmodel$expandNodeNames(paste0("pd[",i,",",g,",1:",J.mark[g],"]"))
#       y.mark.nodes <- Rmodel$expandNodeNames(paste0("y.mark[",i,",",g,",1:",J.mark[g],"]"))
#     }else{
#       pd.nodes <- character(0)
#       y.mark.nodes <- character(0)
#     }
#     if(J.sight[g]>0 & K.sight[g]>0){
#       lam.nodes <- Rmodel$expandNodeNames(paste0("lam[",i,",",g,",1:",J.sight[g],"]"))
#       y.um.nodes <- Rmodel$expandNodeNames(paste0("y.um[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#       y.unk.nodes <- Rmodel$expandNodeNames(paste0("y.unk[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#       lam.um.nodes <- Rmodel$expandNodeNames(paste0("lam.um[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#       lam.unk.nodes <- Rmodel$expandNodeNames(paste0("lam.unk[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#       if(i<=nimbuild$n.marked.all){
#         y.mID.nodes <- Rmodel$expandNodeNames(paste0("y.mID[",i,",",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#         y.mnoID.nodes <- Rmodel$expandNodeNames(paste0("y.mnoID[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#         lam.mnoID.nodes <- Rmodel$expandNodeNames(paste0("lam.mnoID[",g,",1:",J.sight[g],",1:",K.sight[g],"]"))
#       }else{
#         y.mID.nodes <- character(0)
#         y.mnoID.nodes <- character(0)
#         lam.mnoID.nodes <- character(0)
#       }
#     }else{
#       lam.nodes <- character(0)
#       y.mID.nodes <- character(0)
#       y.mnoID.nodes <- character(0)
#       y.um.nodes <- character(0)
#       y.unk.nodes <- character(0)
#       lam.mnoID.nodes <- character(0)
#       lam.um.nodes <- character(0)
#       lam.unk.nodes <- character(0)
#     }
#     conf$addSampler(target=s.target,
#                     type='sSampler3',
#                     control=list(i=i,g=g,xlim=xlim,ylim=ylim,jump.multiplier=2,
#                                  J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,mark.states.all=nimbuild$mark.states,
#                                  n.marked.all=nimbuild$n.marked.all,mark.states=nimbuild$mark.states[i,,],
#                                  s.nodes=s.nodes,pd.nodes=pd.nodes,lam.nodes=lam.nodes,y.mark.nodes=y.mark.nodes,
#                                  y.mID.nodes=y.mID.nodes,y.mnoID.nodes=y.mnoID.nodes,y.um.nodes=y.um.nodes,
#                                  y.unk.nodes=y.unk.nodes,lam.mnoID.nodes=lam.mnoID.nodes,lam.um.nodes=lam.um.nodes,
#                                  lam.unk.nodes=lam.unk.nodes,
#                                  calcNodes=calcNodes),
#                     silent=TRUE)
#   }
# }


#usually a good idea with year-specific sigma
# for(g in 1:n.primary){
#   conf$addSampler(target = c(paste("lam0[",g,"]"),paste("sigma[",g,"]")),
#                   type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
# }

#optional truncated gamma poisson conjugate samplers. 
#I would always use these as long as you keep uniform priors on gamma or gamma[g]
#Typically gives you much greater ESS that propagates to N/N.recruit
#Note: if you add time scaling to model file, need to include that in custom update
#if one gamma per primary occasion
# for(g in 1:(n.primary-1)){
#   target <- paste0("gamma[",g,"]")
#   conf$removeSamplers(target)
#   conf$addSampler(target=target,type=truncGammaPoisSampler,control=list(tau=tau))
# }
#if gamma is fixed
conf$removeSamplers("gamma")
conf$addSampler(target="gamma",type=truncGammaPoisSampler,control=list(tau=tau))

conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc,project=Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 500
plot(mcmc(mvSamples[-c(1:burnin),]))

#reminder what some targets are
data$N
data$N.recruit
data$N.survive
data$N[1] + sum(data$N.recruit) #N.super
sigma.move
rsf.beta 
#check posterior correlations, removing things we can't improve
rem.idx <- c(grep("N",colnames(mvSamples)),
             grep("theta",colnames(mvSamples)))
tmp <- cor(mvSamples[-c(1:burnin),-rem.idx])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)


#Plot N by year with method and mark info
marks.deployed <- colSums(data$mark.deploy) #marks deployed per year
marks.active <- data$n.marked #marks active per year
methods <- ifelse(K.mark > 0 & K.sight > 0, "M-S",
                  ifelse(K.mark > 0, "M",
                         ifelse(K.sight > 0, "S", NA)))

library(vioplot)
vioplot(mvSamples[-c(1:burnin),3:(n.primary+2)],ylim=c(0,200),
        xlim=c(-0.5,n.primary+0.5),ylab="Abundance",line=3)
mtext("Method(s) Used",3,at=0,line=2)
mtext(methods,3,at=1:n.primary,line=2)
mtext("marks deployed",3,at=0,line=1)
mtext(marks.deployed,3,at=1:n.primary,line=1)
mtext("marks active",3,at=0,line=0)
mtext(marks.active,3,at=1:n.primary,line=0)