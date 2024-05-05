NimModel <- nimbleCode({
  D ~ dunif(0,10) #Expected density
  #detection function parameters, shared across sessions
  lam0.fixed ~ dunif(0,10)
  sigma.fixed ~ dunif(0,100)
  
  #sharing genotyping error parameters across sessions
  #genotyping error priors for heterozygote and homozygote loci-level genotypes
  #for each samp.level
  for(sl in 1:samp.levels){
    alpha.het[sl,1:3] <- c(1,1,1)
    p.geno.het[sl,1:3] ~ ddirch(alpha.het[sl,1:3])
    alpha.hom[sl,1:2] <- c(1,1)
    p.geno.hom[sl,1:2] ~ ddirch(alpha.hom[sl,1:2])
  }
  #genotype classification array - shared across sessions
  for(sl in 1:samp.levels){
    theta[1:n.cov,sl,1:max.levels,1:max.levels] <- getTheta(ptype = ptype[1:n.cov,1:max.levels,1:max.levels],
                                                            p.geno.het = p.geno.het[sl,1:3],
                                                            p.geno.hom = p.geno.hom[sl,1:2], n.levels=n.levels[1:n.cov])
  }
  
  for(g in 1:N.session){ #sessions
    lambda[g] <- D*area[g] #expected N
    N[g] ~ dpois(lambda[g]) #realized N
    #put session-specific priors here if not shared
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    
    #genotype frequency priors - these *are* session-specific
    for(m in 1:n.cov){
      for(k in 1:n.levels[m]){
        alpha[g,m,k] <- 1 #dirichlet prior parameters
      }
      gammaMat[g,m,1:n.levels[m]] ~ ddirch(alpha[g,m,1:n.levels[m]])
    }
    
    for(i in 1:M[g]){
      for(m in 1:n.cov){
        G.true[g,i,m]~dcat(gammaMat[g,m,1:n.levels[m]]) #True genotypes
      }
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      y.true[g,i,1:J[g]] ~ dPoissonVector(lam=lam[g,i,1:J[g]]*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod
    }
    
    #genotype observation process, vectorized over reps
    for(l in 1:n.samples[g]){
      for(m in 1:n.cov){
        #custom distribution to skip missing values. Model won't work if you sample unobserved data.
        G.obs[g,l,m,1:n.rep[g]] ~ dcat2(theta=theta[m,samp.type[g,l],G.true[g,ID[g,l],m],1:n.levels[m]],
                                        n.rep = n.rep[g],na.ind=na.ind[g,l,m,1:n.rep[g]]) 
      }
    }
    #calculate number of inds captured
    capcounts[g,1:M[g]] <- Getcapcounts(y.true=y.true[g,1:M[g],1:J[g]]) #intermediate object
    #must use ID and G.latent somewhere to make nimble happy. Sticking them here, not used in function.
    n[g] <- Getncap(capcounts=capcounts[g,1:M[g]],ID=ID[g,1:n.samples[g]],G.latent=G.latent[g,1:M[g],1:n.cov])
  }
})# end model
