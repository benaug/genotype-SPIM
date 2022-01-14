#1 locus

library(snow)
library(doSNOW)
library(foreach)
cores=2
reps=2
cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
registerDoSNOW(cl.tmp)

out=foreach(rep=1:reps) %dopar% {
  library(coda)
  library(nimble)
  source("sim.genoSPIM.poisson.R")
  source("build.genos.R")
  source("init.data.poisson.R")
  source("map.genos.R")
  source("NimbleModel genoSPIM Poisson.R")
  source("Nimble Functions genoSPIM.R")
  
  #make sure to run this line or the MCMC sampler will not work!
  nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
  nimbleOptions('MCMCjointlySamplePredictiveBranches') 
  
  load("unique.genos.RData")
  
  n.levels=unlist(lapply(unique.genos,nrow)) #how many loci-level genotypes per locus?
  
  #This function creates objects to determine which classifications are 1) correct, 2) false allele, and 3) allelic dropout
  built.genos=build.genos(unique.genos)
  ptype=built.genos$ptype #list of length n.cov with each element being an n.levels[m] x n.levels[m] matrix contianing 
  
  #OK, moving on to the genotype frequencies.
  load("gammameans.RData")#loci-level genotype frequency estimates for fisher data set
  
  #Normal SCR stuff
  N=78
  lam0=0.25
  sigma=0.50
  K=10 #number of capture occasoins
  buff=3 #state space buffer. Should be at least 3 sigma.
  X<- expand.grid(3:11,3:11) #trapping array
  n.cov=9 #number of loci, 9 in full data set, but can simulate fewer loci if you want
  n.rep=2 #number of PCR reps per sample.
  
  IDcovs=vector("list",n.cov)
  for(i in 1:n.cov){
    IDcovs[[i]]=1:nrow(unique.genos[[i]])
  }
  gamma=vector("list",n.cov)
  for(i in 1:n.cov){
    # gamma[[i]]=rep(1/n.levels[i],n.levels[i]) #This simulates equal genotype frequencies
    gamma[[i]]=gammameans[[i]] #This uses the frequencies estimated from fisher data set
  }
  pID=rep(1,n.cov) #loci-level sample by replication amplification probabilities (controls level of missing scores)
  p.geno.het=c(0.85,0.149,0.001) #P(correct, allelic dropout,false allele) for heterozygotes (using fisher ests here)
  p.geno.hom=c(0.999,0.001) #P(correct,false allele) for homozygotes
  
  
  data=sim.genoSPIM.poisson(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,#cap-recap parms
                            n.cov=n.cov,pID=pID,n.rep=n.rep,
                            p.geno.hom=p.geno.hom,p.geno.het=p.geno.het,
                            gamma=gamma,IDcovs=IDcovs,ptype=ptype)
  
  #Data augmentation level
  M=150
  J=nrow(X)
  K1D=rep(K,J) #trap operation matrix, number of occasions trap j is operable
  
  #set some gamma inits. Using equal across locus-level genotypes here so we don't use truth
  #note, gamma is a ragged matrix for use in nimble.
  gammaMat=matrix(0,nrow=n.cov,ncol=max(n.levels))
  for(l in 1:n.cov){
    gammaMat[l,1:n.levels[l]]=rep(1/n.levels[l],n.levels[l])
  }
  
  #provide some ballpark inits to initialize latent variables and other data structures
  inits=list(lam0=1,sigma=1,gammaMat=gammaMat,p.geno.het=c(0.95,0.025,0.025),p.geno.hom=c(0.95,0.05)) #plug in some ballpark estimates to initialize data
  nimbuild=init.data.poisson(data=data,M=M,inits=inits)
  
  #book keeping for case where n.rep=1, or n.cov=1
  #nimble does not like structures with any dimensions of 1. Padding the structures
  #in this case lets us use the same code no matter the dimensions.
  if(n.rep==1){
    n.rep.use=2
  }else{
    n.rep.use=n.rep
  }
  if(n.cov==1){
    n.cov.use=2
  }else{
    n.cov.use=n.cov
  }

  #We can't use a list (easily) in nimble, so we use a ragged array instead
  ptypeArray = built.genos$ptypeArray

  #inits for nimble
  Niminits <- list(z=rep(1,M),s=nimbuild$s,G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                   y.true=nimbuild$y.true,G.latent=nimbuild$G.latent,theta=nimbuild$thetaArray,
                   lam0=inits$lam0,sigma=inits$sigma)

  #constants for Nimble
  J=nrow(data$X)
  constants<-list(M=M,J=J,K=data$K,K1D=K1D,n.samples=nimbuild$n.samples,n.cov=n.cov.use,n.rep=n.rep.use,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,na.ind=nimbuild$G.obs.NA.indicator,
                  n.levels=n.levels,max.levels=max(n.levels),ptype=built.genos$ptypeArray)

  #supply data to nimble
  Nimdata<-list(y.true=matrix(NA,nrow=M,ncol=J),G.obs=nimbuild$G.obs,
                G.true=matrix(NA,nrow=M,ncol=n.cov.use),ID=rep(NA,nimbuild$n.samples),
                z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))


  # set parameters to monitor
  parameters<-c('psi','lam0','sigma','N','n','p.geno.het','p.geno.hom')

  # Build the model, configure the mcmc, and compile
  # can ignore warnings about 1) ID in constants 2) possible size mismatch for G.obs.
  start.time<-Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
  conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE)

  #conf$printSamplers() #shows the samplers used for each parameter and latent variable

  ###Two *required* sampler replacements

  ##Here, we remove the default sampler for y.true
  #and replace it with the custom "IDSampler".
  conf$removeSampler("G.obs")
  conf$removeSampler("y.true")
  conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                  type = 'IDSampler',control = list(M=M,J=J,K1D=K1D,n.cov=n.cov.use,n.samples=nimbuild$n.samples,
                                                    n.rep=n.rep.use,this.j=nimbuild$this.j,G.obs=data$G.obs,
                                                    na.ind=nimbuild$G.obs.NA.indicator,n.levels=n.levels),
                  silent = TRUE)

  #replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
  conf$removeSampler("G.true")
  for(i in 1:M){
    for(m in 1:n.cov.use){ #don't need to update first cat bc it is mark status
      conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
                      type = 'GSampler',
                      control = list(i = i,m=m,n.levels=n.levels,n.rep=n.rep.use,
                                     na.ind=nimbuild$G.obs.NA.indicator[,m,]), silent = TRUE)
    }
  }

  # ###Two *optional* sampler replacements:
  #replace default activity center sampler that updates x and y locations separately with a joint update
  conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
  for(i in 1:M){
    # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
    #                 type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""), #do not adapt covariance bc s's not deterministically linked to unmarked individuals
                    type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=250),silent = TRUE)
  }

  #block update for lam0 and sigma helps when data sparse enough to cause correlated posteriors.
  #Slice seems more efficient that RW_block, but slower and may depend on data sparsity
  #If including covariate effects on lam0 or sigma, maybe start with default samplers and see
  #which parameters' posteriors are correlated before deciding what, if anything, to block.
  conf$removeSampler(c("lam0","sigma"))
  conf$addSampler(target = c(paste("lam0"),paste("sigma")),
                  type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)

  # Build and compile
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  # Run the model. Can ignore nimble errors about G.obs value NA or NaN, due to padding to keep dimensions constant for nimble
  start.time2<-Sys.time()
  Cmcmc$run(20000,reset=FALSE) #can extend run by rerunning this line
  end.time<-Sys.time()
  end.time-start.time  # total time for compilation, replacing samplers, and fitting
  end.time-start.time2 # post-compilation run time

  mvSamples = as.matrix(Cmcmc$mvSamples)

  return(mvSamples)
}
stopCluster(cl.tmp)