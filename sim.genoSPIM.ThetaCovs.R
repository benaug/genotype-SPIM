e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#simulates genotyping error rates as function of 1 continuous covariate

sim.genoSPIM.ThetaCovs <-
  function(N=NA,lam0=NA,sigma=NA,theta.d=NA,K=NA,alpha0.het=NA,
           alpha1.het=NA,alpha0.hom=NA,alpha1.hom=NA,
           X=NA,buff=NA,n.cov=NA,n.rep=NA,
           pID=NA,gamma=NA,IDcovs=NA,ptype=NA,obstype="poisson"){
    #error checks
    if(length(gamma)!=n.cov)stop("gamma must be of length n.cov")
    for(l in 1:n.cov){
      if(length(gamma[[l]])!=length(IDcovs[[l]]))stop("gamma[[l]] must have one element per element of IDcovs[[l]]")
      if(sum(gamma[[l]])!=1)stop("gamma[[l]] must sum to 1")
    }
    if(length(alpha0.het)!=2|length(alpha1.het)!=2)warning("alpha0.het and alpha1.het must both be of length 2 unless using SNPs")
    if(length(alpha0.hom)!=1|length(alpha1.hom)!=1)stop("alpha0.hom and alpha1.hom must both be of length 1")
    
    # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    
    #simulate IDcovs
    n.levels=unlist(lapply(gamma,length))
    G.true=matrix(NA,nrow=N,ncol=n.cov) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:n.cov){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    # Capture individuals
    y <-array(0,dim=c(N,J,K))
    if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rpois(1,lamd[i,j])
          }
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(theta.d))stop("Must provide overdispersion parameter theta.d for negbin obsmod")
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rnbinom(1,mu=lamd[i,j],size=theta.d)
          }
        }
      } 
    }else{
      stop("Observation model not recognized.")
    }
    
    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught=which(apply(y,c(1),sum)>0)
    y=y[caught,,]
    G.true=G.true[caught,]
    if(n.cov==1){
      G.true=matrix(G.true,ncol=1)
    }
    s=s[caught,]
    if(K==1){
      y=array(y,dim=c(dim(y),1))
    }
    n=length(caught)
    n.samples=sum(y)
    G.cap=matrix(NA,nrow=n.samples,ncol=n.cov)
    ID=rep(1:nrow(y),times=rowSums(y))
    
    idx=1
    A=array(0,dim=c(dim(y),n.samples))  #configuration matrix: indicator matrix for which individual i occassion j  trap k corresponds to sample l. used to convert corrupt IDcovs to corrupt capture history
    for(i in 1:length(caught)){ #loop through inds (uncaptured already removed)
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y[i,j,k]){ #then samples
              G.cap[idx,]=G.true[i,]
              A[i,j,k,idx]=1
              idx=idx+1
            }
          }
        }
      }
    }
    ycap=aperm(apply(A,c(2,3,4),sum),c(3,1,2))
    
    #double check that observed data can reconstruct true data
    yobs=array(0,dim=c(max(ID),J,K))
    for(i in 1:n.samples){
      map=which(A[,,,i]==1,arr.ind=TRUE)
      yobs[ID[i],map[2],map[3]]=yobs[ID[i],map[2],map[3]]+1
    }
    
    ycheck=array(0,dim=dim(y))
    for(i in 1:n.samples){
      idx2=which(A[,,,i]>0,arr.ind=TRUE)
      ycheck[idx2[1],,]=ycheck[idx2[1],,]+A[idx2[1],,,i]
    }
    if(!all(ycheck==y)){
      stop("Error in data simulator!")
    }

    #build theta with sample covariates
    samp.cov=rnorm(n.samples,0,1)
    useSNPs=FALSE
    if(length(alpha0.het)==2){ #microsats
      logodds.het=matrix(NA,nrow=n.samples,ncol=2)
      p.geno.het=matrix(NA,nrow=n.samples,ncol=3)
      logodds.hom=rep(NA,n.samples)
      p.geno.hom=matrix(NA,nrow=n.samples,ncol=2)
      denom.het=denom.hom=rep(NA,n.samples)
      for(l in 1:n.samples){
        for(i in 1:2){
          logodds.het[l,i] <- alpha0.het[i] + samp.cov[l]*alpha1.het[i]
        }
        denom.het[l] <- 1 + sum(exp(logodds.het[l,1:2]))
        p.geno.het[l,2:3] <- exp(logodds.het[l,1:2])/denom.het[l]
        p.geno.het[l,1] <- 1/denom.het[l]
        logodds.hom[l] <- alpha0.hom + samp.cov[l]*alpha1.hom
        #homozytote multinomial (logistic) regression (2 outcomes)
        denom.hom[l] <- 1 + exp(logodds.hom[l])
        p.geno.hom[l,2] <- exp(logodds.hom[l])/denom.hom[l]
        p.geno.hom[l,1] <- 1/denom.hom[l]
      }
    }else if(length(alpha0.het)==1){#SNPs
      useSNPs=TRUE
      logodds.het=rep(NA,n.samples)
      p.geno.het=matrix(NA,nrow=n.samples,ncol=2)
      logodds.hom=rep(NA,n.samples)
      p.geno.hom=matrix(NA,nrow=n.samples,ncol=2)
      denom.het=denom.hom=rep(NA,n.samples)
      for(l in 1:n.samples){
        logodds.het[l] <- alpha0.het + samp.cov[l]*alpha1.het
        denom.het[l] <- 1 + exp(logodds.het[l])
        p.geno.het[l,1] <- exp(logodds.het[l])/denom.het[l]
        p.geno.het[l,2] <- 1/denom.het[l]
        logodds.hom[l] <- alpha0.hom + samp.cov[l]*alpha1.hom
        denom.hom[l] <- 1 + exp(logodds.hom[l])
        p.geno.hom[l,1] <- exp(logodds.hom[l])/denom.hom[l]
        p.geno.hom[l,2] <- 1/denom.hom[l]
      }
    }

    theta=vector("list",n.cov)
    for(m in 1:n.cov){
      theta[[m]]=array(0,dim=c(n.samples,n.levels[m],n.levels[m]))
      for(l2 in 1:n.samples){
        for(l in 1:n.levels[m]){
          if(!any(ptype[[m]][l,]==2)){#homozygote
            theta[[m]][l2,l,which(ptype[[m]][l,]==1)]=p.geno.hom[l2,1]
            theta[[m]][l2,l,which(ptype[[m]][l,]==3)]=p.geno.hom[l2,2]*(1/sum(ptype[[m]][l,]==3))
          }else{
            theta[[m]][l2,l,which(ptype[[m]][l,]==1)]=p.geno.het[l2,1]
            theta[[m]][l2,l,which(ptype[[m]][l,]==2)]=p.geno.het[l2,2]*(1/sum(ptype[[m]][l,]==2))
            if(!useSNPs){
              theta[[m]][l2,l,which(ptype[[m]][l,]==3)]=p.geno.het[l2,3]*(1/sum(ptype[[m]][l,]==3))
            }
          }
        }
      }
    }
    
    #observation error
    G.error=array(NA,dim=c(n.samples,n.cov,n.rep))
    for(k in 1:n.rep){
      G.error[,,k]=G.cap
      for(l in 1:n.cov){
        for(i in 1:n.samples){
          if(rbinom(1,1,pID[l])==1){
            for(j in 1:n.levels[l]){
              if(G.cap[i,l]==j){
                G.error[i,l,k]=sample(IDcovs[[l]],1,prob=theta[[l]][i,j,])
              }
            }
          }else{
            G.error[i,l,k]=NA
          }
        }
      }
    }
    #find errors that occurred
    G.Obstype=array(0,dim=c(n.samples,n.cov,n.rep))
    for(k in 1:n.rep){
      for(l in 1:n.cov){
        for(i in 1:n.samples){
          if(is.na(G.error[i,l,k]))next#is missing
          if(G.error[i,l,k]==G.cap[i,l]){#is correct
            G.Obstype[i,l,k]=1
          }else{#is an error
            G.Obstype[i,l,k]=ptype[[l]][G.error[i,l,k],G.cap[i,l]]
          }
        }
      }
    }
    getmode = function(v) {
      uniqv = unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    #generate a crude consensus genotype
    G.consensus=apply(G.error,c(1,2),function(x){getmode(x[x!=0])})
    #how many of the crude consensus genotypes are corrupted?
    corrupted=sum(apply(G.consensus==G.cap,1,function(x){any(x==FALSE)}),na.rm=TRUE)
    
    #get this.j, this.k
    tmp=t(apply(ycap,1,function(x){which(x==1,arr.ind=TRUE)}))
    this.j=tmp[,1]
    this.k=tmp[,2]
    
    out=list(y=y,this.j=this.j,this.k=this.k,G.true=G.true,G.obs=G.error,n.cov=n.cov,n.levels=n.levels,
             n.samples=length(this.j),IDlist=list(n.cov=n.cov,IDcovs=IDcovs,ptype=ptype),
             ID=ID,X=X,K=K,buff=buff,s=s,n=nrow(y),corrupted=corrupted,G.Obstype=G.Obstype,obstype=obstype,
             p.geno.het=p.geno.het,p.geno.hom=p.geno.hom,samp.cov=samp.cov)
    
    return(out)
  }