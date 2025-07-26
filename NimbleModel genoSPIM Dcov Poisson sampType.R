NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #detection parameters
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
 
  #Density model
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space
  
  for(i in 1:M){
    for(m in 1:n.cov){
      G.true[i,m] ~ dcat(gammaMat[m,1:n.levels[m]])
    }
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~  dunif(xlim[1],xlim[2])
    s[i,2] ~  dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
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
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],G.latent=G.latent[1:M,1:n.cov])
})# end model
