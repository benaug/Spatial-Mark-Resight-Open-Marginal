library(nimble)
library(coda)
source("sim.JS.SMR.Dcov.Generalized.R")
source("init.SMR.Dcov.Open.Generalized.R")
source("Nimble Model JS-SMR-Dcov Generalized.R")
source("Nimble Functions JS-SMR-Dcov Generalized.R") #contains custom distributions and updates
source("sSampler Dcov Open Marginal Generalized.R") # activity center sampler that proposes from prior when z.super=0.
source("mask.check.R")
#must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

n.primary <- 6 #number of primary occasions
phi <- rep(0.8,n.primary-1) #per-capita recruitment by primary occasion
gamma <- rep(0.2,n.primary-1) #per-capita recruitment by primary occasion
tau <- rep(1,n.primary-1) #duration of each primary-occasion interval
p0 <- rep(0.25,n.primary) #marking process p0
lam0 <- rep(0.25,n.primary) #sighting process lam0
sigma <- rep(0.5,n.primary) #detection function scale by primary occasion
p.mark <- rep(0.75,n.primary) #probability of marking given captured in marking process by primary occasion
#Number of occasions per primary occasion per method
#to skip sampling by a method in a primary occasion, set its K=0
K.mark <- c(5,0,0,5,0,0) #marking occasions by primary occasion
K.sight <- c(5,5,5,5,5,5) #resighting occasions by primary occasion
if(length(K.mark)!=length(K.sight))stop("K.mark and K.sight must be same length")
if(length(K.mark)!=n.primary)stop("K.mark and K.sight must be of length n.primary")

#theta is probability of observing each sample type for marked and unmarked individuals
#assuming the same over primary occasions
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)

#make an SCR trapping array. Making the trapping array size vary by session
#For occasions with no marking or sighting, insert a trap matrix with 0 rows, matrix(0,nrow=0,ncol=2)
X.sight <- vector("list",n.primary)
X.sight[[1]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[2]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[3]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[4]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[5]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[6]] <- as.matrix(expand.grid(1:10,1:10))

X.mark <- vector("list",n.primary)
X.mark[[1]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[2]] <- matrix(0,nrow=0,ncol=2)
X.mark[[3]] <- matrix(0,nrow=0,ncol=2)
X.mark[[4]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[5]] <- matrix(0,nrow=0,ncol=2)
X.mark[[6]] <- matrix(0,nrow=0,ncol=2)

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

### Habitat covariate stuff###
#get x and y extent for each grid separately, then merge
xlim <- ylim <- matrix(NA,n.primary,2)
buff <- 2 #state space buffer around traps
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
set.seed(154)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(50,50),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
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
D.beta1 <- 1
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
data <- sim.JS.SMR.Dcov.Generalized(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,tau=tau,
                                    InSS=InSS,phi=phi,gamma=gamma,n.primary=n.primary,
                                    theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                                    p0=p0,lam0=lam0,sigma=sigma,K.mark=K.mark,K.sight=K.sight,
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

data$N #abundance by primary occasion
colSums(apply(data$y.mark>0,c(1,2),sum)>0) #total marking process captures per primary occasion
colSums(data$mark.deploy) #total marks deployed per primary occasion
rowSums(data$mark.deploy) #total marks deployed per captured individual
data$n.marked #marks active per primary occasion

#total detected individuals
colSums(apply(data$truth$y,c(1,2),sum)>0)
#marked spatial recaps
table(apply(1*(data$truth$y.mark>0),c(1,2),sum))

#visualize all realized activity centers
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X.sight.all,pch=4,lwd=2)
points(X.mark.all,pch=4,col="darkred",lwd=2)
points(data$truth$s,pch=16)

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
  points(data$truth$s[data$truth$z[,g]==1,1],data$truth$s[data$truth$z[,g]==1,2],pch=16) #activity centers
  if(data$n.marked[g]>0){
    for(i in 1:data$n.marked[g]){
      id <- data$ID.marked[[g]][i]
      traps <- matrix(numeric(0),nrow=0,ncol=2)
      if(data$J.sight[g]>0){
        trapcaps <- which(data$y.mID[id,g,]>0)
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
      s <- data$s[id,]
      points(s[1],s[2],col="goldenrod",pch=16)
      if(nrow(traps)>0){
        for(j in 1:nrow(traps)){
          lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
        }
      }
      tel.idx <- which(data$tel.ID == id)
      if(length(tel.idx) > 0){
        tel.g.idx <- which(data$tel.session[tel.idx,]==g)
        if(length(tel.g.idx) > 0){
          nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
          if(nloc>0){
            for(l in 1:nloc){
              lines(x=c(s[1],data$locs[tel.idx,tel.g.idx,l,1]),
                    y=c(s[2],data$locs[tel.idx,tel.g.idx,l,2]),
                    col="gray80")
            }
            points(data$locs[tel.idx,tel.g.idx,1:nloc,1],data$locs[tel.idx,tel.g.idx,1:nloc,2],
                   pch = 16, cex = 0.5, col = "lightblue")
            points(s[1],s[2],col="darkblue",pch = 16)
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

#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
inits <- list(p0=rep(0.1,n.primary),lam0=rep(0.25,n.primary),sigma=rep(0.5,n.primary)) #initializing with 1 parameter per session, just set all to same value
#This function structures the simulated data to fit the model in Nimble (some more restructing below)
nimbuild <- init.SMR.Dcov.Open.Generalized(data,inits,M=M)

#plot to check s inits by primary occasion
for(g in 1:n.primary){
  image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),
        main=paste("primary occasion",g),xlab="X",ylab="Y",col=cols1)
  points(data$X.sight[[g]],pch=4,lwd=2)
  points(data$X.mark[[g]],pch=4,lwd=2,col="darkred")
  points(nimbuild$s[nimbuild$z[,g]==1,1],nimbuild$s[nimbuild$z[,g]==1,2],pch=16) #initialized activity centers
  if(n.marked[g]>0){
    for(i in 1:n.marked[g]){
      id <- nimbuild$ID.marked[i,g]
      traps <- matrix(numeric(0),nrow=0,ncol=2)
      if(J.sight[g]>0){
        trapcaps <- which(nimbuild$y.mID[id,g,1:J.sight[g]]>0)
        if(length(trapcaps)>0){
          traps <- rbind(traps,nimbuild$X.sight[g,trapcaps,1:2])
        }
      }
      if(J.mark[g]>0){
        trapcaps2 <- which(nimbuild$y.mark[id,g,1:J.mark[g]]>0)
        if(length(trapcaps2)>0){
          traps <- rbind(traps,nimbuild$X.mark[g,trapcaps2,1:2])
        }
      }
      s <- nimbuild$s[nimbuild$ID.marked[i,g],]
      points(s[1],s[2],col="goldenrod",pch=16)
      if(nrow(traps)>0){
        for(j in 1:nrow(traps)){
          lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
        }
      }
    }
  }
}

#these indicate in which primary occasion marking/sighting occurs and how many total sessions of each
mark.g <- which(K.mark!=0)
sight.g <- which(K.sight!=0)
n.mark.g <- length(mark.g)
n.sight.g <- length(sight.g)

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.primary=n.primary,M=M,J.mark=J.mark,J.sight=J.sight,xlim=xlim,ylim=ylim,
                  K1D.mark=nimbuild$K1D.mark,K1D.sight=nimbuild$K1D.sight,
                  D.cov=D.cov,cellArea=cellArea,n.cells=n.cells,res=res,
                  n.marked.all=nimbuild$n.marked.all,tau=tau,
                  n.tel.sessions=data$n.tel.sessions,tel.session=data$tel.session,max.n.tel.locs=max.n.tel.locs,
                  tel.ID=data$tel.ID,n.tel.inds=data$n.tel.inds,n.locs.ind=data$n.locs.ind,
                  mark.g=mark.g,sight.g=sight.g,n.mark.g=n.mark.g,
                  n.sight.g=n.sight.g)
#inits for Nimble
Niminits <- list(N=nimbuild$N,N.survive=nimbuild$N.survive,N.recruit=nimbuild$N.recruit,
                 ER=nimbuild$N.recruit,N.super=nimbuild$N.super,z.super=nimbuild$z.super,
                 z=nimbuild$z,z.start=nimbuild$z.start,z.stop=nimbuild$z.stop,
                 s=nimbuild$s,phi.fixed=0.5,D0=nimbuild$N[1]/(sum(InSS)*res^2),D.beta1=0,
                 p0=inits$p0[mark.g],lam0=inits$lam0[sight.g],sigma.fixed=inits$sigma[1])

#data for Nimble
Nimdata <- list(y.mark=nimbuild$y.mark, #marking process
                y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                mark.states=nimbuild$mark.states, #mark state history (who is marked in each primary occasion)
                tel.z.states=nimbuild$tel.z.states, #telemetry z state observations
                dummy.data=nimbuild$dummy.data,cells=cells,InSS=InSS,
                X.mark=nimbuild$X.mark,X.sight=nimbuild$X.sight,locs=data$locs)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','N.super','lambda.y1',
                'phi.fixed','p0','lam0','sigma.fixed','theta.marked','theta.unmarked',
                'D0','D.beta1')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('phi.fixed','gamma','p0','lam0','sigma.fixed','theta.marked','theta.unmarked[2:3]')
#use this above if theta.unmarked is primary occasion specific (no change for theta.marked): paste('theta.unmarked[1:',n.primary,',2:3',']')
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,nodes=config.nodes)

#add N/z sampler
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration?
#20% of M seems reasonable, but optimal will depend on data set
y.mark.nodes <- Rmodel$expandNodeNames(paste0("y.mark[1:",M,",1:",n.primary,",1:",max(J.mark),"]"))
pd.nodes <- Rmodel$expandNodeNames(paste0("pd[1:",M,",1:",n.primary,",1:",max(J.mark),"]"))
lam.nodes <- Rmodel$expandNodeNames(paste0("lam[1:",M,",1:",n.primary,",1:",max(J.sight),"]"))
y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[1:",n.primary,",1:",max(J.sight),"]"))
y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[1:",n.primary,",1:",max(J.sight),"]"))
lam.um.nodes <- Rmodel$expandNodeNames(paste("lam.um[1:",n.primary,",1:",max(J.sight),"]"))
lam.unk.nodes <- Rmodel$expandNodeNames(paste("lam.unk[1:",n.primary,",1:",max(J.sight),"]"))
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.primary-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.primary-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.primary-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
tel.z.states.nodes <- Rmodel$expandNodeNames(paste0("tel.z.states[1:",M,",1]"))
calcNodes <- c(N.nodes,N.recruit.nodes,y.mark.nodes,y.um.nodes,y.unk.nodes,z.nodes,tel.z.states.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.marked.all=n.marked.all,n.cap.all=n.cap.all,
                                                 n.primary=n.primary,J.mark=J.mark,J.sight=J.sight,
                                                 mark.g=mark.g,sight.g=sight.g,
                                                 n.mark.g=n.mark.g,
                                                 n.sight.g=n.sight.g,
                                                 mark.states=nimbuild$mark.states,
                                                 mark.states.all=nimbuild$mark.states,
                                                 tel.z.states=nimbuild$tel.z.states,
                                                 z.super.ups=z.super.ups,y2D=nimbuild$y2D,
                                                 y.mark.nodes=y.mark.nodes,pd.nodes=pd.nodes,
                                                 y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                 lam.nodes=lam.nodes,tel.z.states.nodes=tel.z.states.nodes,
                                                 lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                 N.nodes=N.nodes,z.nodes=z.nodes,ER.nodes=ER.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 calcNodes=calcNodes), silent = TRUE)

#activity center sampler. This sampler tunes activity center MH proposals when z.super[i]=1 and
#draws from the prior otherwise. Also implements Herliansyah et al. (2024) efficiency improvement
#that does not require summation over all individuals' lam values when updating 1 activity center.
for(i in 1:M){
  loc.nodes <- c()
  tel.idx <- which(data$tel.ID==i)
  if(length(tel.idx)>0){
    for(g in 1:data$n.tel.sessions[tel.idx]){
      loc.nodes <- c(loc.nodes,Rmodel$expandNodeNames(paste0("locs[",tel.idx,",",g,",1:",data$n.locs.ind[tel.idx,g],",1:2]")))
    }
  }
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerDcov',control=list(i=i,res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,
                                                     xlim=xlim,ylim=ylim,J.mark=J.mark,J.sight=J.sight,n.marked.all=nimbuild$n.marked.all,
                                                     n.primary=n.primary,sight.g=sight.g,n.sight.g=n.sight.g,loc.nodes=loc.nodes,
                                                     mark.states=nimbuild$mark.states[i,],mark.states.all=nimbuild$mark.states),
                  silent = TRUE)
}

#usually a good idea with primary occasion-specific sigma
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
Cmcmc$run(2000,reset=FALSE) #can extend run by rerunning this line
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

#check posterior correlations, removing things we can't improve
rem.idx <- c(grep("N",colnames(mvSamples)),
             grep("theta",colnames(mvSamples)))
tmp <- cor(mvSamples[-c(1:burnin),-rem.idx])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)


#Plot N by primary occasion with method and mark info
marks.deployed <- colSums(data$mark.deploy) #marks deployed per primary occasion
marks.active <- data$n.marked #marks active per primary occasion
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

