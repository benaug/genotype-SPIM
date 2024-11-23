NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  sigma ~ dunif(0,100)
  #genotype frequency priors
  for(m in 1:n.cov){
    for(k in 1:n.levels[m]){
      alpha[m,k] <- 1 #dirichlet prior parameters
    }
    gammaMat[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]])
  }
  #genotyping error priors for heterozygote and homozygote loci-level genotypes
  #for each samp.level
  for(sl in 1:samp.levels){
    alpha.het[sl,1:3] <- c(1,1,1)
    p.geno.het[sl,1:3] ~ ddirch(alpha.het[sl,1:3])
    alpha.hom[sl,1:2] <- c(1,1)
    p.geno.hom[sl,1:2] ~ ddirch(alpha.hom[sl,1:2])
  }
 
  #--------------------------------------------------------------
  #likelihoods (except for s priors)
  #--------------------------------------------------------------
  for(i in 1:M){
    z[i] ~ dbern(psi)
    for(m in 1:n.cov){
      G.true[i,m] ~ dcat(gammaMat[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
  #genotype classification array
  for(sl in 1:samp.levels){
    theta[1:n.cov,sl,1:max.levels,1:max.levels] <- getTheta(ptype = ptype[1:n.cov,1:max.levels,1:max.levels],
                                                         p.geno.het = p.geno.het[sl,1:3],
                                                         p.geno.hom = p.geno.hom[sl,1:2], n.levels=n.levels[1:n.cov])
  }
  
  #genotype observation process, vectorized over reps
  for(l in 1:n.samples){
    for(m in 1:n.cov){
      #custom distribution to skip missing values. Won't work if you sample unobserved data.
      G.obs[l,m,1:n.rep] ~ dcat2(theta=theta[m,samp.type[l],G.true[ID[l],m],1:n.levels[m]],
                                 n.rep = n.rep,na.ind=na.ind[l,m,1:n.rep]) 
    }
  }
  #calculate number of inds captured and abundance
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J]) #intermediate object
  #must use ID and G.latent somewhere to make nimble happy. Sticking them here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],G.latent=G.latent[1:M,1:n.cov])
  N <- sum(z[1:M])
})# end model
