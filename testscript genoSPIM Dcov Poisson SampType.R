library(coda)
library(nimble)
library(abind)
source("sim.genoSPIM.Dcov.sampType.R")
source("build.genos.R")
source("init.data.Dcov.poisson.sampType.R")
source("map.genos.R")
source("NimbleModel genoSPIM Dcov Poisson sampType.R")
source("Nimble Functions genoSPIM Dcov Poisson sampType.R")
source("sSampler Dcov.R")

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
lam0 <- 0.25
sigma <- 0.50
K <- 10 #number of capture occasions
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11) #trapping array

#more IDcov stuff
n.rep <- 2 #number of PCR reps per sample.
IDcovs <- vector("list",n.cov)
for(i in 1:n.cov){
  IDcovs[[i]] <- 1:nrow(unique.genos[[i]])
}
gamma <- vector("list",n.cov)
for(i in 1:n.cov){
  # gamma[[i]] <- rep(1/n.levels[i],n.levels[i]) #This simulates equal genotype frequencies
  gamma[[i]] <- gammameans[[i]] #This uses the frequencies estimated from fisher data set
}

samp.levels <- 2 #number of sample type covariates. Each type has it's own genotyping error rates.
#sample by replication amplification probabilities (controls level of missing scores)
pID <- c(0.999,0.25) #one for each sample type in this data simulator
p.geno.het <- vector("list",samp.levels)
p.geno.hom <- vector("list",samp.levels)
#P(correct, allelic dropout,false allele) for heterozygotes (using fisher ests here)
p.geno.het[[1]] <- c(0.806,0.185,0.009) #high quality samples
p.geno.het[[2]] <- c(0.489,0.496,0.015) #low quality samples
#P(correct,false allele) for homozygotes
p.geno.hom[[1]] <- c(0.994,0.006) #high quality samples
p.geno.hom[[2]] <- c(0.999,0.001) #low quality samples
p.sampType <- c(0.52,0.48) #from fisher data, 52% high quality, 48% low

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")


### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
# set.seed(13210) #pretty good one
set.seed(13216)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(1000,1000),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]>12] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]>12] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -1
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X,pch=4,cex=0.75)

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(13435341) #change seed for new data set
data <- sim.genoSPIM.Dcov.sampType(D.beta0=D.beta0,D.beta1=D.beta1,
                                         D.cov=D.cov,InSS=InSS,
                                         lam0=lam0,obstype="poisson",
                                         sigma=sigma,K=K,X=X,
                                         xlim=xlim,ylim=ylim,res=res,
                                         n.cov=n.cov,pID=pID,n.rep=n.rep,
                                         p.geno.hom=p.geno.hom,p.geno.het=p.geno.het,
                                         p.sampType=p.sampType,
                                         gamma=gamma,IDcovs=IDcovs,ptype=ptype)

points(data$s.all,pch=16)

#The observed data are 
#1) the trap and occasion of every "count member". E.g., a count of 3 
#(3 genetic samples at a site-occ) has 3 count members
head(data$this.j)
head(data$this.k)
#2) observed genotype replicates for every sample. Can have missing data indicated with NA
t(data$G.obs[1,,]) #observed genotypes for 1st count member

#3) the sample type covariates for every sample
data$samp.type

#Can use the map function to see which genotypes the enumerated genotypes correspond to
ind <- 1 #change ind number ot look at different individuals
data$G.true[ind,] #True genotype of individual 1, enumerated
map.genos(data$G.true[ind,],unique.genos) #converted back to actual genotypes
#can compare observed genotypes to true genotypes above
#allelic dropout events are when heterozygotes are observed as homozygotes. false allele is any other error.
these.samps <- which(data$ID==ind)
for(i in 1:length(these.samps)){
  print(map.genos(t(data$G.obs[these.samps[i],,]),unique.genos))
}

##Structure simulated data for nimble

#Data augmentation level
M <- 150
J <- nrow(X)
K1D <- rep(K,J) #trap operation matrix, number of occasions trap j is operable
data$K1D <- K1D #add to data object to build data below

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

inits <- list(lam0=1,sigma=1,gammaMat=gammaMat,p.geno.het=p.geno.het.init,p.geno.hom=p.geno.hom.init) #plug in some ballpark estimates to initialize data
nimbuild <- init.data.Dcov.poisson.sampType(data=data,M=M,inits=inits)

#book keeping for case where n.rep=1, or n.cov=1
#nimble does not like structures with any dimensions of 1. Padding the structures
#in this case lets us use the same code no matter the dimensions.
if(n.rep==1){
  n.rep.use <- 2
}else{
  n.rep.use <- n.rep
}
if(n.cov==1){
  n.cov.use <- 2
}else{
  n.cov.use <- n.cov
}
n.samples <- data$n.samples

#We can't use a list (easily) in nimble, so we use a ragged array instead
ptypeArray <- built.genos$ptypeArray

#inits for nimble
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z),D0=sum(nimbuild$z)/(sum(data$InSS)*data$res^2),D.beta1=0,
                 s=nimbuild$s,G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,G.latent=nimbuild$G.latent,theta=nimbuild$thetaArray,
                 lam0=inits$lam0,sigma=inits$sigma)

#constants for Nimble
#here, you probably want to center your D.cov. The one I simulated for this testscript is already centered.
# D.cov.use <- data$D.cov - mean(data$D.cov) #plug this into constants$D.cov if centering
J <- nrow(data$X)
constants <- list(M=M,J=J,K1D=K1D,n.samples=n.samples,n.cov=n.cov.use,n.rep=n.rep.use,
                na.ind=nimbuild$G.obs.NA.indicator,samp.levels=samp.levels,
                n.levels=n.levels,max.levels=max(n.levels),ptype=built.genos$ptypeArray,
                D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res)

#supply data to nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(G.obs=nimbuild$G.obs,ID=rep(NA,n.samples),
              X=as.matrix(data$X),samp.type=data$samp.type,
              dummy.data=dummy.data,cells=cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('lambda.N','lam0','sigma','N','n','p.geno.het','p.geno.hom','gammaMat',
                'D0','D.beta1')

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
config.nodes <- c("p.geno.het","p.geno.hom","psi","sigma","lam0","gammaMat",'D0','D.beta1')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,nodes=config.nodes) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
# conf$removeSampler("G.obs") #nimble will assign sampler here if any missing data. Remove it.
# conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K1D=K1D,n.cov=n.cov.use,n.samples=n.samples,
                                                  n.rep=n.rep.use,this.j=nimbuild$this.j,G.obs=data$G.obs,
                                                  na.ind=nimbuild$G.obs.NA.indicator,n.levels=n.levels,
                                                  samp.type=data$samp.type),
                silent = TRUE)

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
# conf$removeSampler("G.true")
# # this is the "safe" version. It will use a lot of RAM, but will be correct if any parameters
# # depend on G.true besides G.obs, which seems unlikely for genotypes. But if, say, G.true[,1] is "sex", and
# # you specify that sigma varies by sex, this update is correct and the more efficient one below will not be.
# for(i in 1:M){
#   for(m in 1:n.cov.use){
#     conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
#                     type = 'GSampler',
#                     control = list(i = i,m=m,n.levels=n.levels,n.rep=n.rep.use,
#                                    na.ind=nimbuild$G.obs.NA.indicator[,m,],
#                                    samp.type=data$samp.type), silent = TRUE)
#   }
# }
#this is the low RAM version. No parameters can depend on G.true except G.obs
#identify G.true nodes here. Must be in matrix with individuals down rows and loci across columns.
#This update only works with "reps" vectorized in bugs code. Must modify this sampler if you unvectorize those.
G.true.nodes <- Rmodel$expandNodeNames(paste0("G.true[1:",M,",1:",n.cov.use,"]"))
G.obs.nodes <- Rmodel$expandNodeNames(paste0("G.obs[1:",n.samples,",1:",n.cov.use,",1:",n.rep.use,"]"))
calcNodes <- c(G.true.nodes,G.obs.nodes)
conf$addSampler(target = paste0("G.true[1:",M,",1:",n.cov,"]"),
                type = 'GSampler2',
                control = list(M=M, n.cov=n.cov.use,n.levels=n.levels,n.rep=n.rep.use,
                               na.ind=nimbuild$G.obs.NA.indicator,n.samples=nimbuild$n.samples,
                               samp.type=data$samp.type,
                               G.true.nodes=G.true.nodes,G.obs.nodes=G.obs.nodes,
                               calcNodes=calcNodes), silent = TRUE)

z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
#nodes used for update
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 y.nodes=y.nodes,lam.nodes=lam.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)


#RW block update for s[i,1:2]
#only tuned for when z=1. When z=0, it draws from the prior, inhomogenous Poisson PP. 
# conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                 xlim=data$xlim,ylim=data$ylim),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#block update for lam0 and sigma helps when data sparse enough to cause correlated posteriors. 
# conf$removeSampler(c("lam0","sigma")) #often better to keep independent samplers with RW_block
conf$addSampler(target = c("lam0","sigma"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
# conf$removeSampler(c("D0","D.beta1")) #AF_slice may be better for D.covs depending on data set
conf$addSampler(target = c("D0","D.beta1"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
#Can ignore nimble warnings about NA or NaN in ptype and theta
#Can ignore nimble warnings about G.obs value NA or NaN, due to padding to keep dimensions constant for nimble
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

#remove gammaMat posteriors (not that interesting and tons of them) and plot
idx <- grep("gammaMat",colnames(mvSamples))
plot(mcmc(mvSamples[200:nrow(mvSamples),-idx]))

data$N
data$n #number of individuals captured to compare to posterior for n. No uncertainty with enough genotype info.

##Explore ID posteriors
#Assuming ID posterior was monitored in mvSamples2.
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
nrow(mvSamples2) #Need enough posterior iterations for reliable inference. If not, reduce thinning and/or run longer
idx <- grep("ID",colnames(mvSamples2))
plot(mcmc(mvSamples2[2:nrow(mvSamples2),idx]))

library(MCMCglmm)
burnin <- 50
IDpost <- round(posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx])))
#For simulated data sets, comparing posterior mode ID to truth.
#Numbers will not be the same, but all samples with same true ID will have
#same ID in posterior mode when posterior mode is exactly correct. Numbers just don't match up.
cbind(data$ID,round(IDpost))

#calculate posterior probability of pairwise sample matches
#P(sample x belongs to same individual as sample y)
burnin <- 50 #where to start. Don't start at 1, is NA.
n.iter <- nrow(mvSamples2)-burnin+1
pair.probs <- matrix(NA,n.samples,n.samples)
for(i in 1:n.samples){
  for(j in 1:n.samples){
    count <- 0
    for(iter in burnin:n.iter){
      count <- count+1*(mvSamples2[iter,idx[j]]==mvSamples2[iter,idx[i]])
    }
    pair.probs[i,j] <- count/(n.iter-burnin+1)
  }
}

this.samp <- 1 #sample number to look at
pair.probs[this.samp,] #probability this sample is from same individual as all other samples
pair.probs[this.samp,data$ID==data$ID[this.samp]] #for simulated data, these are the other samples truly from same individual


#inspect G.true (true genotype) posteriors
idx <- grep("G.true",colnames(mvSamples2))

#posterior mode of true genotypes. Note, this is the posterior mode of each loci individually
#I expect this should usually be the same at the posterior mode complete genotype, but this
#can also be calculated from this posterior
library(MCMCglmm)
burnin <- 50
G.mode <- round(posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx])))
G.mode <- matrix(G.mode,nrow=M)
#rearrange all G.true samples to look at range of values instead of just the mode
G.samps <- mvSamples2[burnin:nrow(mvSamples2),idx]
G.samps <- array(t(G.samps),dim=c(M,n.cov,nrow(G.samps)))

#look at posterior mode genotype of each individual (not numbered the same as in true data,
#but numbers are the same as in IDpost)
ind <- 1 #change ind number to look at different individuals
G.mode[ind,] #True genotype of focal individual, enumerated
map.genos(G.mode[ind,],unique.genos) #converted back to actual genotypes

#which samples were most commonly assigned to this individual? (assumes you calculated IDpost above)
these.samps <- which(IDpost==ind)
if(length(these.samps>0)){
  for(i in 1:length(these.samps)){
    print(map.genos(t(data$G.obs[these.samps[i],,]),unique.genos))
  }
}else{
  "No sample's posterior mode was this individual"
}

#here we can look at the entire posterior of true genotypes for this individual
#Note, individuals with samples strongly linked to them will have precisely
#estimated true genotypes while individuals without samples strongly linked
#will hae very imprecisely estimated true genotypes. If no samples ever allocate,
#you are just drawing true genotypes from the estimated population-level genotype frequencies
out <- t(apply(G.samps[ind,,],2,FUN=map.genos,unique.genos))
head(out,10)

loci <- 1
loci <- loci+1
table(out[,loci])
