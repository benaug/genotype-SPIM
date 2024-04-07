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
source("sim.genoSPIM.sampType.multi.R")
source("sim.genoSPIM.sampType.R")
source("init.data.poisson.sampType.multi.R")
source("init.data.poisson.sampType.R")
source("NimbleModel genoSPIM Poisson sampType multi.R")
source("Nimble Functions genoSPIM Poisson sampType multi.R")
source("sSampler Multisession.R")

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
obstype <- "poisson"
N.session <- 3
lam0 <- rep(0.25,N.session)
sigma <- rep(0.5,N.session)
K <- c(5,6,7)
buff <- rep(3,N.session) #state space buffer. Should be at least 3 sigma.
X <- vector("list",N.session) #one trapping array per session
X[[1]] <- as.matrix(expand.grid(3:11,3:11))
X[[2]] <- as.matrix(expand.grid(3:10,3:10)+1)
X[[3]] <- as.matrix(expand.grid(3:9,3:9)+2)

#get state space areas
area <- get.area(X,buff)
area


D <- 0.4 #expected D
lambda <- D*area #expected N
lambda
N <- rpois(N.session,lambda)#realized N
N

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

data <- sim.genoSPIM.sampType.multi(N.session=N.session,N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,#cap-recap parms
                  obstype="poisson",
                  n.cov=n.cov,pID=pID,n.rep=n.rep,
                  p.geno.hom=p.geno.hom,p.geno.het=p.geno.het,
                  p.sampType=p.sampType,
                  gamma=gamma,IDcovs=IDcovs,ptype=ptype)
#this data is structured exactly as the single session sampler within a session

##Structure simulated data for nimble
#Data augmentation level for each session
M <- c(150,140,130)
#initialize multisession data structures
# nimbuild=init.SCR.Multi(data,M)

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
nimbuild <- init.data.poisson.sampType.multi(data=data,M=M,inits=inits) #initialize latent states

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

Niminits <- list(z=nimbuild$z,s=nimbuild$s,G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=capcounts,
                 y.true=nimbuild$y.true,G.latent=nimbuild$G.latent,theta=thetaArray,
                 lam0.fixed=1,sigma.fixed=1,N=N.init)

#constants for Nimble
constants <- list(N.session=N.session,M=M,J=J,K1D=nimbuild$K1D,n.samples=n.samples,n.cov=n.cov.use,n.rep=n.rep.use,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,na.ind=nimbuild$G.obs.NA.indicator,
                n.levels=n.levels,max.levels=max(n.levels),ptype=built.genos$ptypeArray,
                samp.levels=samp.levels,area=area)

#supply data to nimble
M.max <- max(M)
J.max <- max(J)
n.samples.max <- max(n.samples)
Nimdata <- list(y.true=array(NA,dim=c(N.session,M.max,J.max)),G.obs=nimbuild$G.obs,
              G.true=array(NA,dim=c(N.session,M.max,n.cov.use)),ID=matrix(NA,N.session,n.samples.max),
              z=matrix(NA,N.session,M.max),X=nimbuild$X,capcounts=matrix(NA,N.session,M.max),
              samp.type=nimbuild$samp.type)


# set parameters to monitor
parameters <- c('lam0.fixed','sigma.fixed','N','D','n','p.geno.het','p.geno.hom','gammaMat')

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c('ID',"G.true")
nt <- 1 #thinning rate
nt2 <- 50#thin more

# Build the model, configure the mcmc, and compile
# can ignore warnings about 1) ID in constants 2) possible size mismatch for G.obs.
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#can use "nodes" argument in configureMCMC below to omit y.true and G.true that are replaced below for faster
#configuration
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt, monitors2=parameters2,thin2=nt2,useConjugacy = TRUE)

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
conf$removeSampler("G.obs") #nimble will assign sampler here if any missing data. Remove it.
conf$removeSampler("y.true")
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
conf$removeSampler("G.true")
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

conf$removeSampler("s")
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
#Can ignore nimble warnings about NA or NaN in ptype and theta
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

N #true N
unlist(lapply(data,function(x){x$n})) #number of individuals captured to compare to posterior for n. No uncertainty with enough genotype info.
