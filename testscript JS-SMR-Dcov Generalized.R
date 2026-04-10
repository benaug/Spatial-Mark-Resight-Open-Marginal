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

n.year <- 6 #number of years
phi <- rep(0.8,n.year-1) #yearly per-capita recruitment
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
p0 <- rep(0.1,n.year) #marking process p0
lam0 <- rep(0.25,n.year) #sighting process lam0
sigma <- rep(0.5,n.year) #yearly detection function scale
#Number of occasions per year per method
#to skip sampling by a method in a year, set its K=0
K.mark <- c(5,0,0,5,0,0) #yearly marking occasions
K.sight <- c(5,5,5,5,5,5) #yearly resighting occasions
if(length(K.mark)!=length(K.sight))stop("K.mark and K.sight must be same length")
if(length(K.mark)!=n.year)stop("K.mark and K.sight must be of length n.year")

#theta is probability of observing each sample type for marked and unmarked individuals
#assuming the same over years
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)

#make an SCR trapping array. Making the trapping array size vary by session
#I think the code currently requires marking/sighting traps in all sessions
#even when not used. Will fix that.
X.sight <- vector("list",n.year)
X.sight[[1]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[2]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[3]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[4]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[5]] <- as.matrix(expand.grid(1:10,1:10))
X.sight[[6]] <- as.matrix(expand.grid(1:10,1:10))

X.mark <- vector("list",n.year)
X.mark[[1]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[2]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[3]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[4]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[5]] <- as.matrix(expand.grid(3:8,3:8))
X.mark[[6]] <- as.matrix(expand.grid(3:8,3:8))

### Habitat covariate stuff###
#get x and y extent for each grid separately, then merge
xlim <- ylim <- matrix(NA,n.year,2)
buff <- 2 #state space buffer around traps
X.both <- vector("list",n.year)
for(g in 1:n.year){
  X.both[[g]] <- rbind(X.mark[[g]],X.sight[[g]])
  xlim[g,] <- range(X.both[[g]][,1]) + c(-buff,buff)
  ylim[g,] <- range(X.both[[g]][,2]) + c(-buff,buff)
}
xlim <- c(min(xlim[,1]),max(xlim[,2]))
ylim <- c(min(ylim[,1]),max(ylim[,2]))

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim - x.shift
ylim <- ylim - y.shift
for(g in 1:n.year){
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
for(g in 1:n.year){
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
#what is implied expected year 1 N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected year 1 N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density in Year 1",col=cols1)
points(X.sight.all,pch=4,cex=1,lwd=2)
points(X.mark.all,pch=4,cex=1,lwd=2,col="darkred")

#Mark/Telemetry settings - For now, we assume mark history is known and deaths observed if marked at time of death
#this is a simplified scenario for using telemetry collars for marks
n.tel.locs <- 15 #number of locs per individual
mark.year.pars <- c(2,2,3) #parameters for truncated poisson: c(lambda, lower truncation, upper truncation)
#data simulator requires lower bound be 1 or higher. 1 means it fails before 2nd year
#mark lifetime frequencies for mark.year.pars
table(rtruncpois(10000,lambda=mark.year.pars[1],lower=mark.year.pars[2],upper=mark.year.pars[3]))/10000
#marking protocol: #1) never replace a mark if currently collared on capture 2) always replace
mark.protocol <- 2 

# simulate some data
set.seed(390297) #change seed for new data set
data <- sim.JS.SMR.Dcov.Generalized(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,
            InSS=InSS,phi=phi,gamma=gamma,n.year=n.year,
            theta.marked=theta.marked,theta.unmarked=theta.unmarked,
            p0=p0,lam0=lam0,sigma=sigma,K.mark=K.mark,K.sight=K.sight,
            X.mark=X.mark,X.sight=X.sight,xlim=xlim,ylim=ylim,res=res,
            mark.year.pars=mark.year.pars,mark.protocol=mark.protocol,
            n.tel.locs=n.tel.locs)

#what is observed data? Note data objects have all n.years with all 0 data if no effort for a method
#Could be structured without years with no effort, but that would require more work changing custom
#N/z updates.

# str(data$y.mark) #marking process history: n.marked.all x n.year x J.mark.max
# str(data$y.mID) #marked with ID sighting history: n.marked.all x n.year x J.sight.max
# str(data$y.mnoID) #marked with no ID sighting history: n.year x J.sight.max
# str(data$y.um) #unmarked sighting history: n.year x J.sight.max
# str(data$y.unk) #unknown marked status sighting history: n.year x J.sight.max
# str(data$mark.states) #mark status history: n.marked.all x n.year
# str(data$tel.z.states) #telemetry survival observations: n.marked.all x n.year
# str(data$locs) #telemetry locations: max(n.marked[1:n.year]) x n.year x n.tel.locs x 2

data$N #yearly abundance
colSums(apply(data$y.mark>0,c(1,2),sum)>0) #marks deployed per year 
#(actually this ^ is only true if mark.protocol=2, will over count with 
#mark.protocol=1, need to fix that, but not important for fitting model)
data$n.marked #marks active per year

#total detected individuals
colSums(apply(data$truth$y,c(1,2),sum)>0)
#marked spatial recaps
table(apply(1*(data$truth$y.mark>0),c(1,2),sum))

#visualize all realized activity centers
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X.sight.all,pch=4,lwd=2)
points(X.mark.all,pch=4,col="darkred",lwd=2)
points(data$truth$s,pch=16)

#visualize detections by year. only showing SCR and identified SMR detections
for(g in 1:n.year){
  image(data$x.vals,data$y.vals,matrix(data$D.cov*data$InSS,data$n.cells.x,data$n.cells.y),
        main=paste("Year",g),xlab="X",ylab="Y",col=cols1)
  if(K.sight[g]>0){
    points(data$X.sight[[g]],pch=4,lwd=2)
  }
  if(K.mark[g]>0){
    points(data$X.mark[[g]],pch=4,lwd=2,col="darkred")
  }
  points(data$truth$s[data$truth$z[,g]==1,1],data$truth$s[data$truth$z[,g]==1,2],pch=16) #activity centers
  if(data$n.marked[g]>0){
    for(i in 1:data$n.marked[g]){
      id <- data$ID.marked[[g]][i]
      trapcaps <- which(data$y.mID[id,g,]>0)
      traps <-  data$X.sight[[g]][1:data$J.sight[g],][trapcaps,]
      trapcaps2 <- which(data$y.mark[id,g,]>0)
      traps2 <-  data$X.mark[[g]][1:data$J.mark[g],][trapcaps2,]
      traps <- rbind(traps,traps2)
      s <- data$s[id,]
      points(s[1],s[2],col="goldenrod",pch=16)
      if(nrow(traps)>0){
        for(j in 1:nrow(traps)){
          lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
        }
      }
      tel.idx <- which(data$tel.ID == id)
      if(length(tel.idx) > 0){
        tel.g.idx <- which(data$tel.year[tel.idx,]==g)
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
n.year <- data$n.year
n.marked <- data$n.marked
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
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
inits <- list(p0=rep(0.1,n.year),lam0=rep(0.25,n.year),sigma=rep(0.5,n.year)) #initializing with 1 parameter per session, just set all to same value
#This function structures the simulated data to fit the model in Nimble (some more restructing below)
nimbuild <- init.SMR.Dcov.Open.Generalized(data,inits,M=M)

#plot to check s inits by year
for(g in 1:n.year){
  image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),
        main=paste("Year",g),xlab="X",ylab="Y",col=cols1)
  points(data$X.sight[[g]],pch=4,lwd=2)
  points(data$X.mark[[g]],pch=4,lwd=2,col="darkred")
  points(nimbuild$s[nimbuild$z[,g]==1,1],nimbuild$s[nimbuild$z[,g]==1,2],pch=16) #initialized activity centers
  for(i in 1:n.marked[g]){
    trapcaps <- which(nimbuild$y.mID[nimbuild$ID.marked[i,g],g,]>0)
    traps <-  nimbuild$X.sight[g,1:J.sight[g],][trapcaps,]
    trapcaps2 <- which(nimbuild$y.mark[nimbuild$ID.marked[i,g],g,]>0)
    traps2 <-  nimbuild$X.mark[g,1:J.mark[g],][trapcaps2,]
    traps <- rbind(traps,traps2)
    s <- nimbuild$s[nimbuild$ID.marked[i,g],]
    points(s[1],s[2],col="goldenrod",pch=16)
    if(nrow(traps)>0){
      for(j in 1:nrow(traps)){
        lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
      }
    }
  }
}

#these indicate in which year marking/sighting occurs and how many total sessions of each
mark.years <- which(K.mark!=0)
sight.years <- which(K.sight!=0)
n.mark.years <- length(mark.years)
n.sight.years <- length(sight.years)

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.year=n.year,M=M,J.mark=J.mark,J.sight=J.sight,xlim=xlim,ylim=ylim,
                  K1D.mark=nimbuild$K1D.mark,K1D.sight=nimbuild$K1D.sight,
                  D.cov=D.cov,cellArea=cellArea,n.cells=n.cells,res=res,
                  n.marked.all=nimbuild$n.marked.all,
                  n.tel.years=data$n.tel.years,tel.year=data$tel.year,
                  tel.ID=data$tel.ID,n.tel.inds=data$n.tel.inds,n.locs.ind=data$n.locs.ind,
                  mark.years=mark.years,sight.years=sight.years,n.mark.years=n.mark.years,
                  n.sight.years=n.sight.years)
#inits for Nimble
Niminits <- list(N=nimbuild$N,N.survive=nimbuild$N.survive,N.recruit=nimbuild$N.recruit,
                 ER=nimbuild$N.recruit,N.super=nimbuild$N.super,z.super=nimbuild$z.super,
                 z=nimbuild$z,z.start=nimbuild$z.start,z.stop=nimbuild$z.stop,
                 s=nimbuild$s,phi.fixed=0.5,D0=nimbuild$N[1]/(sum(InSS)*res^2),D.beta1=0,
                 p0=inits$p0[mark.years],lam0=inits$lam0[sight.years],sigma.fixed=inits$sigma[1])

#data for Nimble
Nimdata <- list(y.mark=nimbuild$y.mark, #marking process
                y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                mark.states=nimbuild$mark.states, #mark state history (who is marked in each year)
                tel.z.states=nimbuild$tel.z.states, #telemetry z state observations
                dummy.data=nimbuild$dummy.data,cells=cells,InSS=InSS,
                X.mark=nimbuild$X.mark,X.sight=nimbuild$X.sight,locs=data$locs)

# set parameters to monitor
parameters <- c('N','gamma.fixed','N.recruit','N.survive','N.super','lambda.y1',
                'phi.fixed','p0','lam0','sigma.fixed','theta.marked','theta.unmarked',
                'D0','D.beta1')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('phi.fixed','gamma.fixed','p0','lam0','sigma.fixed','theta.marked','theta.unmarked[2:3]')
#use this above if theta.unmarked is year specific (no change for theta.marked): paste('theta.unmarked[1:',n.year,',2:3',']')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = FALSE)

#add N/z sampler
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration?
#20% of M seems reasonable, but optimal will depend on data set
y.mark.nodes <- Rmodel$expandNodeNames(paste0("y.mark[1:",M,",1:",n.year,",1:",max(J.mark),"]"))
pd.nodes <- Rmodel$expandNodeNames(paste0("pd[1:",M,",1:",n.year,",1:",max(J.sight),"]"))
lam.nodes <- Rmodel$expandNodeNames(paste0("lam[1:",M,",1:",n.year,",1:",max(J.sight),"]"))
y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[1:",n.year,",1:",max(J.sight),"]"))
y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[1:",n.year,",1:",max(J.sight),"]"))
lam.mnoID.nodes <- Rmodel$expandNodeNames(paste("lam.mnoID[1:",n.year,",1:",max(J.sight),"]"))
lam.um.nodes <- Rmodel$expandNodeNames(paste("lam.um[1:",n.year,",1:",max(J.sight),"]"))
lam.unk.nodes <- Rmodel$expandNodeNames(paste("lam.unk[1:",n.year,",1:",max(J.sight),"]"))
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
tel.z.states.nodes <- Rmodel$expandNodeNames(paste0("tel.z.states[1:",M,",1]"))
calcNodes <- c(N.nodes,N.recruit.nodes,y.mark.nodes,y.um.nodes,y.unk.nodes,z.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.marked.all=nimbuild$n.marked.all,
                                                 n.year=n.year,J.mark=J.mark,J.sight=J.sight,
                                                 mark.years=mark.years,sight.years=sight.years,
                                                 n.mark.years=n.mark.years,
                                                 n.sight.years=n.sight.years,
                                                 mark.states=nimbuild$mark.states,
                                                 tel.z.states=nimbuild$tel.z.states,
                                                 z.super.ups=z.super.ups,y2D=nimbuild$y2D,
                                                 y.mark.nodes=y.mark.nodes,pd.nodes=pd.nodes,
                                                 y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                 lam.nodes=lam.nodes,lam.mnoID.nodes=lam.mnoID.nodes,
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
    for(g in 1:data$n.tel.years[tel.idx]){
      loc.nodes <- c(loc.nodes,Rmodel$expandNodeNames(paste0("locs[",tel.idx,",",g,",1:",data$n.locs.ind[tel.idx,g],",1:2]")))
    }
  }
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerDcov',control=list(i=i,res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,
                                                     xlim=xlim,ylim=ylim,J.mark=J.mark,J.sight=J.sight,n.marked.all=nimbuild$n.marked.all,
                                                     n.year=n.year,loc.nodes=loc.nodes,mark.states=nimbuild$mark.states[i,]),
                  silent = TRUE)
}

#usually a good idea with year-specific sigma
# for(g in 1:n.year){
#   conf$addSampler(target = c(paste("lam0[",g,"]"),paste("sigma[",g,"]")),
#                   type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
# }

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
plot(mcmc(mvSamples[-c(1:250),]))

#reminder what some targets are
data$N
data$N.recruit
data$N.survive
data$N[1] + sum(data$N.recruit) #N.super

#check posterior correlations, removing things we can't improve
rem.idx <- c(grep("N",colnames(mvSamples)),
             grep("theta",colnames(mvSamples)))
tmp <- cor(mvSamples[-c(1:500),-rem.idx])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)


#Plot N by year with method and mark info
marks.deployed <- colSums(apply(data$y.mark>0,c(1,2),sum)>0) #marks deployed per year
marks.active <- data$n.marked #marks active per year
methods <- ifelse(K.mark > 0 & K.sight > 0, "M-S",
                  ifelse(K.mark > 0, "M",
                         ifelse(K.sight > 0, "S", NA)))

library(vioplot)
vioplot(mvSamples[-c(1:500),3:(n.year+2)],ylim=c(0,200),
        xlim=c(-0.5,n.year+0.5),ylab="Abundance",line=3)
mtext("Method(s) Used",3,at=0,line=2)
mtext(methods,3,at=1:n.year,line=2)
mtext("marks deployed",3,at=0,line=1)
mtext(marks.deployed,3,at=1:n.year,line=1)
mtext("marks active",3,at=0,line=0)
mtext(marks.active,3,at=1:n.year,line=0)
