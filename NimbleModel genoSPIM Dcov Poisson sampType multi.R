NimModel <- nimbleCode({
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
    #Density covariates
    # D.beta0[g] ~ dnorm(0,sd=10)
    D0[g] ~ dunif(0,100) #Density intercept prior. uninformative, diffuse dnorm on log scale can cause neg bias
    D.beta1[g] ~ dnorm(0,sd=10)
    #Density model
    D.intercept[g] <- D0[g]*cellArea[g]
    # D.intercept <- exp(D.beta0)*cellArea
    lambda.cell[g,1:n.cells[g]] <- InSS[g,1:n.cells[g]]*exp(D.beta1[g]*D.cov[g,1:n.cells[g]])
    pi.denom[g] <- sum(lambda.cell[g,1:n.cells[g]])
    pi.cell[g,1:n.cells[g]] <- lambda.cell[g,1:n.cells[g]]/pi.denom[g] #expected proportion of total N in cell c
    #Abundance model
    lambda.N[g] <- D.intercept[g]*pi.denom[g] #Expected N
    N[g] ~ dpois(lambda.N[g]) #realized N in state space
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
        G.true[g,i,m] ~ dcat(gammaMat[g,m,1:n.levels[m]]) #True genotypes
      }
      #dunif() here implies uniform distribution within a grid cell
      #also tells nimble s's are in continuous space, not discrete
      s[g,i,1] ~  dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~  dunif(ylim[g,1],ylim[g,2])
      #get cell s_i lives in using look-up table
      s.cell[g,i] <- cells[g,trunc(s[g,i,1]/res[g])+1,trunc(s[g,i,2]/res[g])+1]
      #categorical likelihood for this cell, equivalent to zero's trick
      #also disallowing s's in non-habitat
      dummy.data[g,i] ~ dCell(pi.cell[g,s.cell[g,i]])
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
    #calculate number of inds captured and abundance
    capcounts[g,1:M[g]] <- Getcapcounts(ID=ID[g,1:n.samples[g]],M=M[g]) #intermediate object
    #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
    n[g] <- Getncap(capcounts=capcounts[g,1:M[g]],G.latent=G.latent[g,1:M[g],1:n.cov])
  }
})# end model
