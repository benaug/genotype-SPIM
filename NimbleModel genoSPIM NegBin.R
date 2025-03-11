NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  theta.d ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  sigma ~ dunif(0,100)
  #genotype frequency priors
  for(m in 1:n.cov){
    for(k in 1:n.levels[m]){
      alpha[m,k] <- 1 #dirichlet prior parameters
    }
    gammaMat[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]])
  }
  #genotyping error priors for heterozygote and homozygote loci-level genotypes
  alpha.het[1:3] <- c(1,1,1)
  p.geno.het[1:3] ~ ddirch(alpha.het[1:3])
  alpha.hom[1:2] <- c(1,1)
  p.geno.hom[1:2] ~ ddirch(alpha.hom[1:2])
 
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
    p[i,1:J] <- theta.d/(theta.d + lam[i,1:J])
    y.true[i,1:J] ~ dNBVector(p=p[i,1:J],theta=theta.d*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of NB RVs is NB with theta=N*theta
  }
  #genotype classification array
  theta[1:n.cov,1:max.levels,1:max.levels] <- getTheta(ptype = ptype[1:n.cov,1:max.levels,1:max.levels],
                                                       p.geno.het = p.geno.het[1:3],
                                                       p.geno.hom = p.geno.hom[1:2], n.levels=n.levels[1:n.cov])
  #genotype observation process, vectorized over reps
  for(l in 1:n.samples){
    for(m in 1:n.cov){
      #custom distribution to skip missing values. Won't work if you sample unobserved data.
      G.obs[l,m,1:n.rep] ~ dcat2(theta=theta[m,G.true[ID[l],m],1:n.levels[m]],
                                 n.rep = n.rep,na.ind=na.ind[l,m,1:n.rep]) 
    }
  }
  #calculate number of inds captured and abundance
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],G.latent=G.latent[1:M,1:n.cov])
  N <- sum(z[1:M])
})# end model
