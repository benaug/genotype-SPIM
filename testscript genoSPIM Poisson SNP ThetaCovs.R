##This version works for SNPs. Set up for fixed genotyping error rates across SNPs (runs much faster).
#will need to modify a lot to allow SNP-level error rates
#This version also allows genotyping error rates to vary by sample covariates.
#This test script is set up to use 1 continuous covaraite, e.g., percent amplification
#per sample

library(coda)
library(nimble)
library(abind)
source("sim.genoSPIM.ThetaCovs.R")
source("build.genos.R")
source("init.data.poisson.R")
source("map.genos.R")
source("NimbleModel genoSPIM Poisson SNP ThetaCovs.R")
source("Nimble Functions genoSPIM Poisson SNP ThetaCovs.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#First, let's structure some SNPs.
n.cov <- 20 #number of SNP loci
unique.genos <- vector("list")
for(m in 1:n.cov){
  unique.genos[[m]] <- matrix(c(1,1,2,2,1,2),nrow=3,byrow=TRUE)
}
n.levels <- rep(3,n.cov) #each SNP has 3 possible genotypes, 11, 12, 22

#This function creates objects to determine which classifications are 1) correct, 2) false allele, and 3) allelic dropout
built.genos <- build.genos(unique.genos)
ptype <- built.genos$ptype #list of length n.cov with each element being an n.levels[m] x n.levels[m] matrix containing 
#indicator for error type for true genotypes along rows and classified genotype along columns


#Normal SCR stuff
N <- 50
lam0 <- 0.25
sigma <- 0.50
K <- 10 #number of capture occasions
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11) #trapping array
n.rep <- 3 #number of PCR reps per sample.

IDcovs <- vector("list",n.cov)
for(i in 1:n.cov){
  IDcovs[[i]] <- 1:nrow(unique.genos[[i]])
}
gamma=vector("list",n.cov)
for(i in 1:n.cov){
  gamma[[i]] <- rep(1/n.levels[i],n.levels[i]) #This simulates equal genotype frequencies
}

#sample by replication amplification probabilities (controls level of missing scores)
pID <- rep(0.95,n.cov) #one for each covariate in this data simulator


#heterozygotes
alpha0.het <- 4 #intercept for logodds(allelic dropout/correct)
alpha1.het <-  3 #coef for logodds(allelic dropout/correct)
#homozygotes
alpha0.hom <- 4 #intercept for logodds(false allele/correct)
alpha1.hom <- 1.5 #coef for logodds(false allele/correct)

data <- sim.genoSPIM.ThetaCovs(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,#cap-recap parms
                              obstype="poisson",
                              n.cov=n.cov,pID=pID,n.rep=n.rep,
                              alpha0.het=alpha0.het,alpha1.het=alpha1.het,
                              alpha0.hom=alpha0.hom,alpha1.hom=alpha1.hom,
                              gamma=gamma,IDcovs=IDcovs,ptype=ptype)

#Look at covariate-error relationships:
par(mfrow=c(1,1),ask=FALSE)
plot(data$p.geno.het[,1]~data$samp.cov,ylim=c(0,1),pch=16,col="darkgreen",main="Heterozygous Loci-level Genotypes",
     xlab="some continouous sample covariate",ylab="P(classification type)") #P(correct classification|heterozygote)
points(data$p.geno.het[,2]~data$samp.cov,pch=16,col="goldenrod") #P(allelic dropout|heterozygote)
legend(x="bottomright",legend=c("Correct","Allelic Dropout"),col=c("darkgreen","goldenrod"),
       pch=c(16,16,16),cex=0.75)
plot(data$p.geno.hom[,1]~data$samp.cov,ylim=c(0,1),pch=16,col="darkgreen",main="Homozygous Loci-level Genotypes",
     xlab="some continouous sample covariate",ylab="P(classification type)") #P(correct classification|heterozygote)
points(data$p.geno.hom[,2]~data$samp.cov,pch=16,col="darkred") #P(false allele|heterozygote)
legend(x="bottomright",legend=c("Correct","False Allele"),col=c("darkgreen","darkred","goldenrod"),
       pch=c(16,16,16),cex=0.75)



#The observed data are 
#1) the trap and occasion of every "count member". E.g., a count of 3 
#(3 genetic samples at a site-occ) has 3 count members
head(data$this.j)
head(data$this.k)
#2) observed genotype replicates for every sample. Can have missing data indicated with NA
t(data$G.obs[1,,]) #observed genotypes for 1st count member

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
M <- 125
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

#Theta Covs note: for simplicity, we can use the same data initializer and just 
#initialize theta to be the same regardless of sample covariates
inits <- list(lam0=1,sigma=1,gammaMat=gammaMat,p.geno.het=c(0.95,0.05),p.geno.hom=c(0.95,0.05)) #plug in some ballpark estimates to initialize data
nimbuild <- init.data.poisson(data=data,M=M,inits=inits)


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

#Just pull out the first one since they're all the same for SNPs
ptypeMatrix  <-  built.genos$ptypeArray[1,,]
# thetaMatrix <- nimbuild$thetaArray[1,,,] #do i need this with covariates?

#inits for nimble
Niminits <- list(z=rep(1,M),s=nimbuild$s,G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,G.latent=nimbuild$G.latent, #theta=thetaMatrix,
                 lam0=inits$lam0,sigma=inits$sigma, #using lam0 and sigma truth for inits. dont do in practice.
                 alpha0.het=alpha0.het,alpha1.het=alpha1.het,alpha0.hom=alpha0.hom,alpha1.hom=alpha1.hom) 
#setting alpha inits to true values. Don't do in practice, but might need ballpark inits 
#if very little replication or you could get label switching.


#constants for Nimble
J <- nrow(data$X)
constants <- list(M=M,J=J,K1D=K1D,n.samples=n.samples,n.cov=n.cov.use,n.rep=n.rep.use,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,na.ind=nimbuild$G.obs.NA.indicator,
                ptype=ptypeMatrix,samp.cov=data$samp.cov)

#supply data to nimble
Nimdata <- list(y.true=matrix(NA,nrow=M,ncol=J),G.obs=nimbuild$G.obs,
              G.true=matrix(NA,nrow=M,ncol=n.cov.use),ID=rep(NA,n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))


# set parameters to monitor
parameters <- c('psi','lam0','sigma','N','n','alpha0.het','alpha0.hom','alpha1.het','alpha1.hom')

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c('ID',"G.true",'gammaMat') #monitoring the IDs and the true genotypes so we can explore them below
nt <- 1 #thinning rate
nt2 <- 50 #thin more

# Build the model, configure the mcmc, and compile
# can ignore warnings about 1) ID in constants 2) possible size mismatch for G.obs.
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#using config.nodes for faster compilation-skip nimble assigning samplers that need to be removed and replaced
#if you add parameters to the model, add to config.nodes to assign them samplers
config.nodes <- c('alpha0.het','alpha0.hom','alpha1.het','alpha1.hom',"psi","sigma","lam0","gammaMat","z")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,
                      nodes=config.nodes,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
# conf$removeSampler("G.obs") #nimble will assign sampler here if any missing data. Remove it.
# conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K1D=K1D,n.cov=n.cov.use,n.samples=n.samples,
                                                  n.rep=n.rep.use,this.j=nimbuild$this.j,G.obs=data$G.obs,
                                                  na.ind=nimbuild$G.obs.NA.indicator),
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
#                     control = list(i = i,m=m,n.rep=n.rep.use,samp.type=data$samp.type,
#                                    na.ind=nimbuild$G.obs.NA.indicator[,m,]), silent = TRUE)
#   }
# }
#this is the low RAM version. No parameters can depend on G.true except G.obs
#identify G.true nodes here. Must be in matrix with individuals down rows and loci across columns.
#This update only works with "reps" vectorized in bugs code. Must modify this sampler if you unvectorize those.
G.true.nodes <- Rmodel$expandNodeNames(paste0("G.true[1:",M,",1:",n.cov.use,"]"))
G.obs.nodes <- Rmodel$expandNodeNames(paste0("G.obs[1:",n.samples,",1:",n.cov.use,",1:",n.rep.use,"]"))
calcNodes <- c(G.true.nodes,G.obs.nodes)
conf$addSampler(target = paste0("G.true[1:",M,",1:",n.cov.use,"]"),
                type = 'GSampler2',
                control = list(M=M,n.cov=n.cov.use,n.rep=n.rep.use,
                               na.ind=nimbuild$G.obs.NA.indicator,n.samples=nimbuild$n.samples,
                               G.true.nodes=G.true.nodes,G.obs.nodes=G.obs.nodes,samp.type=data$samp.type,
                               calcNodes=calcNodes), silent = TRUE)


#RW block update for s[i,1:2]
#only tuned for when z=1. When z=0, it draws from the prior, assumed to be uniform. 
# conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#block update for lam0 and sigma helps when data sparse enough to cause correlated posteriors. 
# conf$removeSampler(c("lam0","sigma")) #often better to keep independent samplers
conf$addSampler(target = c("lam0","sigma"),
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
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

plot(mcmc(mvSamples[200:nrow(mvSamples),],))

data$n #number of individuals captured to compare to posterior for n. No uncertainty with enough genotype info.

##Explore ID posteriors
#Assuming ID posterior was monitored in mvSamples2
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
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
#will have very imprecisely estimated true genotypes. If no samples ever allocate,
#you are just drawing true genotypes from the estimated population-level genotype frequencies
out <- t(apply(G.samps[ind,,],2,FUN=map.genos,unique.genos))
head(out,10)

loci <- 1
loci <- loci+1
table(out[,loci])

