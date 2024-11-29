#This multisession sampler is hard coded such that all possible
#genotypes at each locus, and the number of loci don't vary
#across sessions. Changing this requires changing custom updates.
#number of replications can vary by session.
#genotype error parameters currently shared across sessions,
#but genotype frequencies do vary across sessions.
#I haven't tried this version with only 1 locus, but that seems
#unlikely in practice. Code imported from single session may work in this case.

library(coda)
library(nimble)
library(abind)
source("build.genos.R")
source("map.genos.R")
source("mask.check.R")
source("sim.genoSPIM.Dcov.sampType.multi.R")
source("sim.genoSPIM.Dcov.sampType.R")
source("init.data.Dcov.poisson.sampType.multi.R")
source("init.data.Dcov.poisson.sampType.R")
source("NimbleModel genoSPIM Dcov Poisson sampType multi.R")
source("Nimble Functions genoSPIM Dcov Poisson sampType multi.R")
source("sSampler Dcov Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#We will use the fisher data set parameter estimates from Augustine et al. (2020)
#to specify the genetic portion of the simulation settings

#This is a list of locus-level genotypes. Each list element corresponds to a locus
#and contains a matrix of all possible combinations of locus-level genotypes, one per row.
#the order of alleles in a genotype does not matter
load("unique.genos.RData")
str(unique.genos)
#A real data set should enumerate the observed loci-level genotypes in the order listed here,
#e.g. 152-152 is 1, 152-154 is 2, 152-156 is 3, etc.
unique.genos[[1]] #1st loci

n.levels <- unlist(lapply(unique.genos,nrow)) #how many loci-level genotypes per locus?

#This function creates objects to determine which classifications are 1) correct, 2) false allele, and 3) allelic dropout
built.genos <- build.genos(unique.genos)
ptype <- built.genos$ptype #list of length n.cov with each element being an n.levels[m] x n.levels[m] matrix contianing 
#indicator for error type for true genotypes along rows and classified genotype along columns
#note, "ptype" in list form is used in the data simulator, but ptypeArray is a ragged array used in nimble

#####example 1: true homozygote:
#the first genotype at locus 1 is 152-152. Because it is homozygous, there cannot be allelic dropout.
unique.genos[[1]][1,]
#Now, looking at ptype for true locus-level genotype 1 at locus 1:
#we see only a correct classification or a false allele are possible.
ptype[[1]][1,]

#there is only 1 correct classification event, indicated with a 1:
unique.genos[[1]][ptype[[1]][1,]==1,]

#here are the classifications that are the result of false alleles (everything that is not correct)
unique.genos[[1]][ptype[[1]][1,]==3,]

#####example 2: true heterozygote:
#the second genotype at locus 1 is 152-154. Because it is heterozygous, there can be allelic dropout or false alleles
unique.genos[[1]][2,]
#Now, looking at ptype for true locus-level genotype 2 at locus 1:
#we see correct classification, allelic dropout, or false allele are possible.
ptype[[1]][2,]

#there is only 1 correct classification event, indicated with a 1:
unique.genos[[1]][ptype[[1]][2,]==1,]

#here are the classifications that are the result of allelic dropout (152 or 154 can drop out)
unique.genos[[1]][ptype[[1]][2,]==2,]

#here are the classifications that are the result of false alleles
unique.genos[[1]][ptype[[1]][2,]==3,]

#OK, moving on to the genotype frequencies.
load("gammameans.RData")#loci-level genotype frequency estimates for fisher data set
#These are the frequencies of each locus-level genotype at each locus
str(gammameans)

#Now we have the information required to simulate a data set similar to the fisher data set

#First, let's decide how many loci to use
n.cov <- 9
#discard unused information if you don't use them all
if(n.cov!=9){
  for(i in 9:(n.cov+1)){
    gammameans[[i]] <- NULL
    unique.genos[[i]] <- NULL
    ptype[[i]] <- NULL
  }
}
n.levels <- unlist(lapply(unique.genos,nrow)) #update n.levels in case some loci discarded
#now all lists of length "n.cov"
str(gammameans)
str(unique.genos)
str(ptype)
n.levels

#Normal SCR stuff
N.session <- 3
obstype <- "poisson"
lam0 <- rep(0.25,N.session)
sigma <- rep(0.5,N.session)
K <- c(5,6,7)
buff <- rep(3,N.session) #state space buffer. Should be at least 3 sigma.
#make an SCR trapping array. Making the trapping array size vary by session
X <- vector("list",N.session)
X[[1]] <- as.matrix(expand.grid(1:9,1:9))
X[[2]] <- as.matrix(expand.grid(1:8,1:8))
X[[3]] <- as.matrix(expand.grid(1:10,1:10))

#genotype stuff
n.rep <- c(2,2,3) #number of PCR reps per sample in each session
IDcovs <- vector("list",n.cov) #fixed across sessions. all possible loci-level genotypes the same across sessoins
for(i in 1:n.cov){
  IDcovs[[i]] <- 1:nrow(unique.genos[[i]])
}
#gamma... OK, let's put it in a list of lists. Sessions, then loci. 
#But I am plugging in the same values across sessions here.
gamma <- vector("list",N.session)
for(g in 1:N.session){
  gamma[[g]] <- vector("list",n.cov)
  for(i in 1:n.cov){
    gamma[[g]][[i]] <- gammameans[[i]] #This uses the frequencies estimated from fisher data set
  }
}

samp.levels <- 2 #number of sample type covariates. Each type has it's own genotyping error rates.
#sample by replication amplification probabilities (controls level of missing scores)
pID <- c(0.999,0.25) #one for each sample type in this data simulator

#session-specific genotyping error, again list of lists. samp.levels same across sessions. Simulating shared parms across sessions.
p.geno.het <- vector("list",N.session)
p.geno.hom <- vector("list",N.session)
for(g in 1:N.session){
  p.geno.het[[g]] <- vector("list",samp.levels)
  p.geno.hom[[g]] <- vector("list",samp.levels)
  #P(correct, allelic dropout,false allele) for heterozygotes (using fisher ests here)
  p.geno.het[[g]][[1]] <- c(0.806,0.185,0.009) #high quality samples
  p.geno.het[[g]][[2]] <- c(0.489,0.496,0.015) #low quality samples
  #P(correct,false allele) for homozygotes
  p.geno.hom[[g]][[1]] <- c(0.994,0.006) #high quality samples
  p.geno.hom[[g]][[2]] <- c(0.999,0.001) #low quality samples
}

p.sampType <- c(0.52,0.48) #from fisher data, 52% high quality, 48% low


#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- ylim <- matrix(NA,N.session,2)
for(g in 1:N.session){
  xlim[g,] <- range(X[[g]][,1]) + c(-buff[g],buff[g])
  ylim[g,] <- range(X[[g]][,2]) + c(-buff[g],buff[g])
}

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
for(g in 1:N.session){
  x.shift <- xlim[g,1]
  y.shift <- ylim[g,1]
  xlim[g,] <- xlim[g,] - x.shift
  ylim[g,] <- ylim[g,] - y.shift
  X[[g]][,1] <- X[[g]][,1]- x.shift
  X[[g]][,2] <- X[[g]][,2]- y.shift
}

res <- rep(0.25,N.session) #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- y.vals <- dSS <- cells <- vector("list",N.session)
n.cells <- n.cells.x <- n.cells.y <- rep(NA,N.session)
for(g in 1:N.session){
  x.vals[[g]] <- seq(xlim[g,1]+res[g]/2,xlim[g,2]-res[g]/2,res[g]) #x cell centroids
  y.vals[[g]] <- seq(ylim[g,1]+res[g]/2,ylim[g,2]-res[g]/2,res[g]) #y cell centroids
  dSS[[g]] <- as.matrix(cbind(expand.grid(x.vals[[g]],y.vals[[g]])))
  cells[[g]] <- matrix(1:nrow(dSS[[g]]),nrow=length(x.vals[[g]]),ncol=length(y.vals[[g]]))
  n.cells[g] <- nrow(dSS[[g]])
  n.cells.x[g] <- length(x.vals[[g]])
  n.cells.y[g] <- length(y.vals[[g]])
}

#create a density covariate - one for each session
library(geoR)
D.cov <- vector("list",N.session)
#need a simulated landscape with individuals living around traps to be captured
#these are pretty good
D.seeds <- c(13216,13216,13218)
for(g in 1:N.session){
  set.seed(D.seeds[g])
  D.cov.tmp <- grf(n.cells[g],grid=dSS[[g]],cov.pars=c(1000,1000),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
  D.cov.tmp <- as.numeric(scale(D.cov.tmp)) #scale
  par(mfrow=c(1,1),ask=FALSE)
  D.cov[[g]] <- D.cov.tmp
  image(x.vals[[g]],y.vals[[g]],matrix(D.cov[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," D.cov"),xlab="X",ylab="Y",col=cols1)
}

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
InSS <- vector("list",N.session)
for(g in 1:N.session){
  dSS.tmp <- dSS[[g]] - res[g]/2 #convert back to grid locs
  InSS[[g]] <- rep(1,length(D.cov[[g]]))
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  image(x.vals[[g]],y.vals[[g]],matrix(InSS[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Habitat"))
}

#Density covariates
D.beta0 <- rep(-1,N.session)
D.beta1 <- rep(0.5,N.session)
#what is implied expected N in state space?
for(g in 1:N.session){
  lambda.cell <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  print(sum(lambda.cell)) #expected N in state space
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Expected Density"))
  points(X[[g]],pch=4,cex=0.75) #SCR traps
}

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(13435341) #change seed for new data set
data <- sim.genoSPIM.Dcov.sampType.multi(N.session=N.session,D.beta0=D.beta0,D.beta1=D.beta1,
                                   D.cov=D.cov,InSS=InSS,
                                   lam0=lam0,obstype="poisson",
                                   sigma=sigma,K=K,X=X,
                                   xlim=xlim,ylim=ylim,res=res,
                                   n.cov=n.cov,pID=pID,n.rep=n.rep,
                                   p.geno.hom=p.geno.hom,p.geno.het=p.geno.het,
                                   p.sampType=p.sampType,
                                   gamma=gamma,IDcovs=IDcovs,ptype=ptype)

#simulated N per session
unlist(lapply(data,function(x){x$N}))
#SCR-detected individuals by session. 
unlist(lapply(data,function(x){x$n}))

#Visualize activity centers
for(g in 1:N.session){
  lambda.cell <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Expected Density"))
  points(X[[g]],pch=4,cex=0.75)
  points(data[[g]]$s.all,pch=16)
}

for(g in 1:N.session){
  #function to test for errors in mask set up. 
  mask.check(dSS=data[[g]]$dSS,cells=data[[g]]$cells,n.cells=data[[g]]$n.cells,n.cells.x=data[[g]]$n.cells.x,
             n.cells.y=data[[g]]$n.cells.y,res=data[[g]]$res,xlim=data[[g]]$xlim,ylim=data[[g]]$ylim,
             x.vals=data[[g]]$x.vals,y.vals=data[[g]]$y.vals)
}

##Structure simulated data for nimble
#Data augmentation level
M <- c(175,175,200) #one per session

J <- unlist(lapply(X,nrow))
K <- unlist(lapply(data,function(x){x$K}))
for(g in 1:N.session){
  K1D <- rep(K[g],J[g]) #trap operation matrix, number of occasions trap j is operable
  data[[g]]$K1D <- K1D #add to data object to build data below
}

#for these genotyping inits, I am using the same for each session.
#set some gamma inits. Using equal across locus-level genotypes here so we don't use truth
#note, gamma is a ragged matrix for use in nimble.
gammaMat <- matrix(0,nrow=n.cov,ncol=max(n.levels))
for(l in 1:n.cov){
  gammaMat[l,1:n.levels[l]] <- rep(1/n.levels[l],n.levels[l])
}

#provide some ballpark inits to initialize latent variables and other data structures
#to play nice with nimble, we'll store the genotyping error rates in a ragged matrix
p.geno.het.init <- matrix(NA,nrow=2,ncol=3)
p.geno.hom.init <- matrix(NA,nrow=2,ncol=2)
p.geno.het.init[1,] <- c(0.95,0.025,0.025)
p.geno.het.init[2,] <- c(0.95,0.025,0.025)
p.geno.hom.init[1,] <- c(0.95,0.05)
p.geno.hom.init[2,] <- c(0.95,0.05)

#plug in some ballpark estimates to initialize latent states
inits <- list(lam0=rep(1,N.session),sigma=rep(1,N.session),gammaMat=gammaMat,p.geno.het=p.geno.het.init,p.geno.hom=p.geno.hom.init)
nimbuild <- init.data.Dcov.poisson.sampType.multi(data=data,M=M,inits=inits) #initialize latent states

#book keeping for case where n.rep=1, or n.cov=1
#nimble does not like structures with any dimensions of 1. Padding the structures
#in this case lets us use the same code no matter the dimensions.
n.rep.use <- rep(NA,N.session)
for(g in 1:N.session){
  if(n.rep[g]==1){
    n.rep.use[g] <- 2
  }else{
    n.rep.use[g] <- n.rep[g]
  }
}
if(n.cov==1){
  n.cov.use <- 2
}else{
  n.cov.use <- n.cov
}
n.samples <- nimbuild$n.samples

#We can't use a list (easily) in nimble, so we use a ragged array instead
ptypeArray <- built.genos$ptypeArray

#inits for nimble
capcounts <- apply(nimbuild$y.true,c(1,2),sum,na.rm=TRUE)

#data initializer spits out one thetaArray per session, but we only need 1 if the
#genotyping error rates are shared across sessions as I am assuming here
str(nimbuild$thetaArray)
thetaArray <- nimbuild$thetaArray[1,,,,]#pull out one session to use
N.init <- rowSums(nimbuild$z,na.rm=TRUE) #N.init must be consistent with z.init!
D0.init <- N.init/(unlist(lapply(data,function(x){sum(x$InSS)}))*nimbuild$cellArea)

Niminits <- list(z=nimbuild$z,s=nimbuild$s,N=N.init,D0=D0.init,
                 G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=capcounts,
                 y.true=nimbuild$y.true,G.latent=nimbuild$G.latent,theta=thetaArray,
                 lam0.fixed=1,sigma.fixed=1)

#constants for Nimble
constants <- list(N.session=N.session,M=M,J=J,K1D=nimbuild$K1D,n.samples=n.samples,n.cov=n.cov.use,n.rep=n.rep.use,
                na.ind=nimbuild$G.obs.NA.indicator,samp.levels=samp.levels,
                n.levels=n.levels,max.levels=max(n.levels),ptype=built.genos$ptypeArray,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                D.cov=nimbuild$D.cov,cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells,
                res=nimbuild$res)

#supply data to nimble
M.max <- max(M)
J.max <- max(J)
n.samples.max <- max(n.samples)
Nimdata <- list(G.obs=nimbuild$G.obs,ID=matrix(NA,N.session,n.samples.max),
              z=matrix(NA,N.session,M.max),X=nimbuild$X,capcounts=matrix(NA,N.session,M.max),
              samp.type=nimbuild$samp.type,dummy.data=nimbuild$dummy.data,cells=nimbuild$cells,InSS=nimbuild$InSS)


# set parameters to monitor
parameters <- c('lam0.fixed','sigma.fixed','N','D0','D.beta1','n','p.geno.het','p.geno.hom','gammaMat')

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c('ID',"G.true")
nt <- 1 #thinning rate
nt2 <- 50 #thin more

# Build the model, configure the mcmc, and compile
# can ignore warnings about 1) ID in constants 2) possible size mismatch for G.obs.
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#using config.nodes for faster compilation-skip nimble assigning samplers that need to be removed and replaced
#if you add parameters to the model, add to config.nodes to assign them samplers
config.nodes <- c("p.geno.het","p.geno.hom","D0","D.beta1","sigma.fixed","lam0.fixed","gammaMat")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,
                      nodes=config.nodes,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###*required* sampler replacements
#z/N update
z.ups <- round(M*0.25) # how many z proposals per iteration per session?
conf$removeSampler("N")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y.true[",g,",","1:",M[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,y.nodes,lam.nodes)
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}


##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
# conf$removeSampler("G.obs") #nimble will assign sampler here if any missing data. Remove it.
# conf$removeSampler("y.true")
for(g in 1:N.session){
  conf$addSampler(target = paste0("y.true[1:",g,",1:",M[g],",1:",J[g],"]"),
                  type = 'IDSampler',control = list(M=M[g],J=J[g],K1D=nimbuild$K1D[g,1:J[g]],n.cov=n.cov.use,n.samples=n.samples[g],
                                                    n.rep=n.rep.use[g],this.j=nimbuild$this.j[g,1:n.samples[g]],
                                                    G.obs=nimbuild$G.obs[g,1:n.samples[g],1:n.cov,1:n.rep.use[g]],
                                                    na.ind=nimbuild$G.obs.NA.indicator[g,1:n.samples[g],1:n.cov,1:n.rep.use[g]],
                                                    n.levels=n.levels,g=g,
                                                    samp.type=nimbuild$samp.type[g,1:n.samples[g]]),
                  silent = TRUE)
}

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
# conf$removeSampler("G.true")
# # this is the "safe" version. It will use a lot of RAM, but will be correct if any parameters
# # depend on G.true besides G.obs, which seems unlikely for genotypes. But if, say, G.true[,1] is "sex", and
# # you specify that sigma varies by sex, this update is correct and the more efficient one below will not be.
# for(g in 1:N.session){
#   for(i in 1:M[g]){
#     for(m in 1:n.cov.use){
#       conf$addSampler(target = paste("G.true[",g,",",i,",",m,"]", sep=""),
#                       type = 'GSampler',
#                       control = list(g=g,i=i,m=m,n.samples=n.samples[g],n.levels=n.levels,n.rep=n.rep.use[g],
#                                      G.obs=nimbuild$G.obs[g,1:n.samples[g],,1:n.rep[g]],
#                                      na.ind=nimbuild$G.obs.NA.indicator[g,1:n.samples[g],m,],
#                                      samp.type=nimbuild$samp.type[g,1:n.samples[g]]), silent = TRUE)
#     }
#   }
# }
# this is the low RAM version. No parameters can depend on G.true except G.obs
#identify G.true nodes here. Must be in matrix with individuals down rows and loci across columns.
#This update only works with "reps" vectorized in bugs code. Must modify this sampler if you unvectorize those.
for(g in 1:N.session){
  G.true.nodes <- Rmodel$expandNodeNames(paste0("G.true[",g,",1:",M[g],",1:",n.cov.use,"]"))
  G.obs.nodes <- Rmodel$expandNodeNames(paste0("G.obs[",g,",1:",n.samples[g],",1:",n.cov.use,",1:",n.rep.use[g],"]"))
  calcNodes <- c(G.true.nodes,G.obs.nodes)
  conf$addSampler(target = paste0("G.true[",g,",1:",M[g],",1:",n.cov.use,"]"),
                  type = 'GSampler2',
                  control = list(M =M[g], n.cov=n.cov.use,n.levels=n.levels,n.rep=n.rep.use[g],
                                 na.ind=nimbuild$G.obs.NA.indicator[g,1:nimbuild$n.samples[g],,1:n.rep.use[g]],
                                 n.samples=nimbuild$n.samples[g],
                                 samp.type=nimbuild$samp.type[g,1:nimbuild$n.samples[g]],g=g,
                                 G.true.nodes=G.true.nodes,G.obs.nodes=G.obs.nodes,
                                 calcNodes=calcNodes), silent = TRUE)
  
  
}

#RW block update for s[g,i,1:2]
#only tuned for when z=1. When z=0, it draws from the prior, assumed to be inhomogenous Poisson PP. 
# conf$removeSampler("s")
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,J=J[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                   n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#block update for lam0 and sigma helps when data sparse enough to cause correlated posteriors. 
# conf$removeSampler(c("lam0.fixed","sigma.fixed")) #often better to keep independent samplers
conf$addSampler(target = c("lam0.fixed","sigma.fixed"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

#often better to block these 2 sets of parameters
conf$removeSampler(target = c("D0","D.beta1")) #keep independent if using RW_block, remove if AF_slice
for(g in 1:N.session){
  #AF_slice mixes better, but runs more slowly, with reduced runtime a function of how
  #costly it is to evaluate the likelihoods involved.
  #AF_slice slower speed may be worth it for for D covs
  conf$addSampler(target = c(paste0("D0[",g,"]"),paste0("D.beta1[",g,"]")),
                  type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
#Can ignore nimble warnings about NA or NaN in ptype, theta, and G.true
#Can ignore nimble warnings about G.obs value NA or NaN, due to padding to keep dimensions constant for nimble
start.time2 <- Sys.time()
Cmcmc$run(1500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

#remove gammaMat posteriors (not that interesting and tons of them) and plot
idx <- grep("gammaMat",colnames(mvSamples))
plot(mcmc(mvSamples[250:nrow(mvSamples),-idx]))

#truth
unlist(lapply(data,function(x){x$N})) #realized N
unlist(lapply(data,function(x){x$lambda})) #expected N
unlist(lapply(data,function(x){x$n})) #number of individuals captured to compare to posterior for n. No uncertainty with enough genotype info.
