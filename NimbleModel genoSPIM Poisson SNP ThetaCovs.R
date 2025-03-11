NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  sigma ~ dunif(0,100)
  #genotype frequency priors
  for(m in 1:n.cov){
    for(k in 1:3){
      alpha[m,k] <- 1 #dirichlet prior parameters
    }
    gammaMat[m,1:3] ~ ddirch(alpha[m,1:3])
  }
  #genotyping error as function of covariates
  alpha0.het ~ dlogis(0,1) #heterozygote intercept
  alpha1.het ~ dnorm(0,sd=10) #heterozygote slope
  alpha0.hom ~ dlogis(0,1) #homozygote intercept
  alpha1.hom ~ dnorm(0,sd=10) #homozygote slope
  #--------------------------------------------------------------
  #likelihoods (except for s priors)
  #--------------------------------------------------------------
  for(i in 1:M){
    z[i] ~ dbern(psi)
    for(m in 1:n.cov){
      G.true[i,m] ~ dcat(gammaMat[m,1:3])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
  
  #genotype observation process, vectorized over reps
  for(l in 1:n.samples){
    #heterozygote multinomial link (2 outcomes)
    logodds.het[l] <- alpha0.het + samp.cov[l]*alpha1.het
    denom.het[l] <- 1 + exp(logodds.het[l])
    p.geno.het[l,1] <- exp(logodds.het[l])/denom.het[l]
    p.geno.het[l,2] <- 1/denom.het[l]
    logodds.hom[l] <- alpha0.hom + samp.cov[l]*alpha1.hom
    #homozygote multinomial (logistic) link (2 outcomes)
    denom.hom[l] <- 1 + exp(logodds.hom[l])
    p.geno.hom[l,1] <- exp(logodds.hom[l])/denom.hom[l]
    p.geno.hom[l,2] <- 1/denom.hom[l]
    
    #sample-level genotype classification array
    theta[l,1:3,1:3] <- getTheta(ptype = ptype[1:3,1:3],p.geno.het = p.geno.het[l,1:2],p.geno.hom = p.geno.hom[l,1:2])
    for(m in 1:n.cov){
      #custom distribution to skip missing values. Won't work if you sample unobserved data.
      G.obs[l,m,1:n.rep] ~ dcat2(theta=theta[l,G.true[ID[l],m],1:3],
                                 n.rep = n.rep,na.ind=na.ind[l,m,1:n.rep]) 
    }
  }
  #calculate number of inds captured and abundance
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],G.latent=G.latent[1:M,1:n.cov])
  N <- sum(z[1:M])
})# end model
