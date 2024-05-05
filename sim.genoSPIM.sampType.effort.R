e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.genoSPIM.sampType.effort <-
  function(N=NA,beta0.lam0=NA,beta1.lam0=NA,effort=NA,K2D=NA,sigma=NA,theta.d=NA,K=NA,p.geno.het=NA,
           p.geno.hom=NA,X=NA,buff=NA,n.cov=NA,n.rep=NA,p.sampType=NA,
           pID=NA,gamma=NA,IDcovs=NA,ptype=NA,obstype="poisson"){
    #error checks
    if(length(gamma)!=n.cov)stop("gamma must be of length n.cov")
    for(l in 1:n.cov){
      if(length(gamma[[l]])!=length(IDcovs[[l]]))stop("gamma[[l]] must have one element per element of IDcovs[[l]]")
      if(sum(gamma[[l]])!=1)stop("gamma[[l]] must sum to 1")
    } 
    samp.levels=length(p.sampType)
    for(i in 1:samp.levels){
      if(sum(p.geno.hom[[i]])!=1)stop("p.geno.hom must sum to 1 for each sampType")
      if(sum(p.geno.het[[i]])!=1)stop("p.geno.het must sum to 1 for each sampType")
      if(length(p.geno.hom[[i]])!=2)stop("p.geno.hom must be of length 2 for each sampType")
      if(length(p.geno.het[[i]])!=3)warning("p.geno.het must be of length 3 unless simulating SNPs")
    }
    if(length(pID)!=samp.levels)stop("pID must be of length samp.levels")
    
    # simulate a population of activity centers
    s <- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D <- e2dist(s,X)
    J <- nrow(X)
    lam0 <- exp(beta0.lam0 + beta1.lam0*effort)
    lamd <- array(0,dim=c(N,J,K))
    for(j in 1:J){
      for(k in 1:K){
        lamd[,j,k] <- lam0[j,k]*exp(-D[,j]^2/(2*sigma*sigma))
      }
    }
    
    #simulate IDcovs
    n.levels <- unlist(lapply(gamma,length))
    G.true <- matrix(NA,nrow=N,ncol=n.cov) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:n.cov){
        G.true[i,j] <- sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    # Capture individuals
    y <- array(0,dim=c(N,J,K))
    if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k] <- rpois(1,lamd[i,j,k]*K2D[j,k])
          }
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(theta.d))stop("Must provide overdispersion parameter theta.d for negbin obsmod")
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k] <- rnbinom(1,mu=lamd[i,j,k],size=theta.d*K2D[j,k])
          }
        }
      } 
    }else{
      stop("Observation model not recognized.")
    }
    
    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught <- which(apply(y,c(1),sum)>0)
    y <- y[caught,,]
    G.true <- G.true[caught,]
    if(n.cov==1){
      G.true <- matrix(G.true,ncol=1)
    }
    s <- s[caught,]
    if(K==1){
      y <- array(y,dim=c(dim(y),1))
    }
    n <- length(caught)
    n.samples <- sum(y)
    G.cap <- matrix(NA,nrow=n.samples,ncol=n.cov)
    ID <- rep(1:nrow(y),times=rowSums(y))
    
    idx <- 1
    A <- array(0,dim=c(dim(y),n.samples))  #configuration matrix: indicator matrix for which individual i occassion j  trap k corresponds to sample l. used to convert corrupt IDcovs to corrupt capture history
    for(i in 1:length(caught)){ #loop through inds (uncaptured already removed)
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y[i,j,k]){ #then samples
              G.cap[idx,] <- G.true[i,]
              A[i,j,k,idx] <- 1
              idx <- idx+1
            }
          }
        }
      }
    }
    ycap <- aperm(apply(A,c(2,3,4),sum),c(3,1,2))
    
    #double check that observed data can reconstruct true data
    yobs <- array(0,dim=c(max(ID),J,K))
    for(i in 1:n.samples){
      map <- which(A[,,,i]==1,arr.ind=TRUE)
      yobs[ID[i],map[2],map[3]] <- yobs[ID[i],map[2],map[3]]+1
    }
    
    ycheck <- array(0,dim=dim(y))
    for(i in 1:n.samples){
      idx2 <- which(A[,,,i]>0,arr.ind=TRUE)
      ycheck[idx2[1],,] <- ycheck[idx2[1],,]+A[idx2[1],,,i]
    }
    if(!all(ycheck==y)){
      stop("Error in data simulator!")
    }
    #get sample covariates
    samp.type <- sample(samp.levels,nrow(G.cap),replace=TRUE,prob=p.sampType)
    
    theta <- vector("list",n.cov)
    for(m in 1:n.cov){
      theta[[m]] <- array(0,dim=c(samp.levels,n.levels[m],n.levels[m]))
      for(sl in 1:samp.levels){
        for(l in 1:n.levels[m]){
          if(!any(ptype[[m]][l,]==2)){#homozygote
            theta[[m]][sl,l,which(ptype[[m]][l,]==1)] <- (p.geno.hom[[sl]][1])
            theta[[m]][sl,l,which(ptype[[m]][l,]==3)] <- (p.geno.hom[[sl]][2])*(1/sum(ptype[[m]][l,]==3))
          }else{
            theta[[m]][sl,l,which(ptype[[m]][l,]==1)] <- (p.geno.het[[sl]][1])
            theta[[m]][sl,l,which(ptype[[m]][l,]==2)] <- (p.geno.het[[sl]][2])*(1/sum(ptype[[m]][l,]==2))
            theta[[m]][sl,l,which(ptype[[m]][l,]==3)] <- (p.geno.het[[sl]][3])*(1/sum(ptype[[m]][l,]==3))
          }
        }
      }
    }
    
    #observation error
    G.error <- array(NA,dim=c(n.samples,n.cov,n.rep))
    for(k in 1:n.rep){
      G.error[,,k] <- G.cap
      for(l in 1:n.cov){
        for(i in 1:n.samples){
          if(rbinom(1,1,pID[samp.type[i]])==1){
            for(j in 1:n.levels[l]){
              if(G.cap[i,l]==j){
                G.error[i,l,k] <- sample(IDcovs[[l]],1,prob=theta[[l]][samp.type[i],j,])
              }
            }
          }else{
            G.error[i,l,k] <- NA
          }
        }
      }
    }
    #find errors that occurred
    G.Obstype <- array(0,dim=c(n.samples,n.cov,n.rep))
    for(k in 1:n.rep){
      for(l in 1:n.cov){
        for(i in 1:n.samples){
          if(is.na(G.error[i,l,k]))next#is missing
          if(G.error[i,l,k]==G.cap[i,l]){#is correct
            G.Obstype[i,l,k] <- 1
          }else{#is an error
            G.Obstype[i,l,k] <- ptype[[l]][G.error[i,l,k],G.cap[i,l]]
          }
        }
      }
    }
    getmode = function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    #generate a crude consensus genotype
    G.consensus <- apply(G.error,c(1,2),function(x){getmode(x[x!=0])})
    #how many of the crude consensus genotypes are corrupted?
    corrupted <- sum(apply(G.consensus==G.cap,1,function(x){any(x==FALSE)}),na.rm=TRUE)
    
    #get this.j, this.k
    tmp <- t(apply(ycap,1,function(x){which(x==1,arr.ind=TRUE)}))
    this.j <- tmp[,1]
    this.k <- tmp[,2]
    
    out <- list(y=y,this.j=this.j,this.k=this.k,G.true=G.true,G.obs=G.error,n.cov=n.cov,n.levels=n.levels,n.rep=n.rep,
             n.samples=length(this.j),IDlist=list(n.cov=n.cov,IDcovs=IDcovs,ptype=ptype),
             ID=ID,X=X,K=K,buff=buff,s=s,n=nrow(y),corrupted=corrupted,G.Obstype=G.Obstype,obstype=obstype,
             samp.type=samp.type,effort=effort,K2D=K2D)
    
    return(out)
  }