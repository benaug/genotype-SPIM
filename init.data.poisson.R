e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.data.poisson=function(data=NA,M=NA,inits=inits){
  library(abind)
  this.j <- data$this.j
  this.k <- data$this.k
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- data$K
  K1D<-data$K1D
  n.cov=data$n.cov
  n.levels=data$n.levels
  IDcovs=data$IDlist$IDcovs
  
  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  n.samples=length(this.j)
  G.obs=data$G.obs
  if(!is.array(G.obs))stop("G.obs must be an array")
  n.rep=dim(G.obs)[3]
  n.cov=dim(G.obs)[2]
  
  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0<- inits$lam0
  sigma<- inits$sigma
  gammaMat=inits$gammaMat
  p.geno.het=inits$p.geno.het
  p.geno.hom=inits$p.geno.hom
  
  if(length(p.geno.het)!=3)warning("p.geno.het init must be of length 3 unless using SNPs")
  if(length(p.geno.hom)!=2)stop("p.geno.hom init must be of length 2")
  if(sum(p.geno.het)!=1)stop("p.geno.het init must sum to 1")
  if(sum(p.geno.hom)!=1)stop("p.geno.hom init must sum to 1")
  
  #initialize G.obs.true, the true sample-level full categorical identities
  #initializing to the most commonly observed sample by category values across observers
  G.obs.true=matrix(0,nrow=n.samples,ncol=n.cov)
  for(i in 1:n.samples){
    for(l in 1:n.cov){
      vec=G.obs[i,l,]
      if(any(is.na(vec))){
        vec=vec[-which(is.na(vec))]
      }
      if(length(vec)>0){
        choose=as.integer(names(which(table(vec)==max(table(vec)))))
        if(length(choose)>1){
          pick=sample(1:length(choose),1)
          choose=choose[pick]
        }
        G.obs.true[i,l]=choose
      }
    }
  }
  
  #restructure observed data for multinomial evaluation. fill in below
  G.obs.mn=vector("list",n.cov)
  for(l in 1:n.cov){
    G.obs.mn[[l]]=array(0,dim=c(n.samples,n.rep,n.levels[l]))
  }
  ptype=data$IDlist$ptype
  #construct thetaArray
  max.levels=max(n.levels)
  thetaArray=array(NA,dim=c(n.cov,max.levels,max.levels))
  for(m in 1:n.cov){
    for(l in 1:n.levels[m]){
      if(!any(ptype[[m]][l,]==2)){#homozygote
        thetaArray[m,l,which(ptype[[m]][l,]==1)]=(p.geno.hom[1])
        thetaArray[m,l,which(ptype[[m]][l,]==3)]=(p.geno.hom[2])*(1/sum(ptype[[m]][l,]==3))
      }else{
        thetaArray[m,l,which(ptype[[m]][l,]==1)]=(p.geno.het[1])
        thetaArray[m,l,which(ptype[[m]][l,]==2)]=(p.geno.het[2])*(1/sum(ptype[[m]][l,]==2))
        thetaArray[m,l,which(ptype[[m]][l,]==3)]=(p.geno.het[3])*(1/sum(ptype[[m]][l,]==3))
      }
    }
  }
  
  #make G constraints for initializing 
  constraints=constraints.init=matrix(1,nrow=n.samples,ncol=n.samples)
  for(i in 1:n.samples){
    for(j in 1:n.samples){
      guys1=which(G.obs.true[i,]!=0)
      guys2=which(G.obs.true[j,]!=0)
      comp=guys1[which(guys1%in%guys2)]
      if(any(G.obs.true[i,comp]!=G.obs.true[j,comp])){
        constraints.init[i,j]=0
      }
    }
  }
  
  #Build y.true
  y.obs=array(0,dim=c(length(this.j),J,K))
  for(l in 1:n.samples){
    y.obs[l,this.j[l],this.k[l]]=1
  }
  y.true=array(0,dim=c(M,J,K))
  y.true2D=matrix(0,nrow=M,ncol=J)
  y.obs2D=apply(y.obs,c(1,2),sum)
  ID=rep(NA,n.samples)
  idx=1
  for(i in 1:n.samples){
    if(idx>M){
      stop("Need to raise M to initialize y.true")
    }
    traps=which(y.obs2D[i,]>0)
    # y.true2D=apply(y.true,c(1,2),sum)
    if(length(traps)==1){
      cand=which(y.true2D[,traps]>0)#guys caught at same traps
    }else{
      cand=which(rowSums(y.true2D[,traps])>0)#guys caught at same traps
    }
    if(length(cand)>0){
      if(length(cand)>1){#if more than 1 ID to match to, choose first one
        cand=cand[1]
      }
      #Check constraint matrix
      cands=which(ID%in%cand)#everyone assigned this ID
      if(all(constraints.init[i,cands]==1)){#focal consistent with all partials already assigned
        y.true[cand,,]=y.true[cand,,]+y.obs[i,,]
        ID[i]=cand
      }else{#focal not consistent
        y.true[idx,,]=y.obs[i,,]
        y.true2D[idx,traps]=y.true2D[idx,traps]+1
        ID[i]=idx
        idx=idx+1
      }
    }else{#no assigned samples at this trap
      y.true[idx,,]=y.obs[i,,]
      y.true2D[idx,traps]=y.true2D[idx,traps]+1
      ID[i]=idx
      idx=idx+1
    }
  }
  
  #Check assignment consistency with constraints
  checkID=unique(ID)
  for(i in 1:length(checkID)){
    idx=which(ID==checkID[i])
    if(!all(constraints[idx,idx]==1)){
      stop("ID initialized improperly")
    }
  }
  y.true2D=apply(y.true,c(1,2),sum)
  known.vector=c(rep(1,max(ID)),rep(0,M-max(ID)))
  
  #Initialize z
  z=1*(apply(y.true2D,1,sum)>0)
  add=M*(0.5-sum(z)/M)
  if(add>0){
    z[sample(which(z==0),add)]=1 #switch some uncaptured z's to 1.
  }
  
  #Optimize starting locations given where they are trapped.
  s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y.true)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  
  #collapse data to 2D
  y.obs=apply(y.obs,c(1,2),sum)
  y.true2D=apply(y.true,c(1,2),sum)
  D=e2dist(s, X)
  lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))
  ll.y=matrix(NA,M,J)
  for(j in 1:J){
    if(K1D[j]==0&sum(y.true2D[,j])>0)stop("K1D inconsistent with capture data")
    ll.y[,j]=dpois(y.true2D[,j],lamd[,j]*z*K1D[j],log=TRUE)
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting likelihood not finite. Try raising sigma and/or lam0")

  #Initialize G.true
  G.true=matrix(0, nrow=M,ncol=n.cov)
  for(i in 1:max(ID)){
    idx=which(ID==i)
    if(length(idx)==1){
      G.true[i,]=G.obs.true[idx,]
    }else{
      if(ncol(G.obs.true)>1){
        G.true[i,]=apply(G.obs.true[idx,],2, max) #consensus
      }else{
        G.true[i,]=max(G.obs.true[idx,])
      }
    }
  }
  #Fill in missing values, create indicator for them
  G.latent=G.true==0
  for(j in 1:n.cov){
    fix=G.true[,j]==0
    G.true[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gammaMat[j,1:n.levels[j]])
  }
  
  this.j=apply(y.obs,1,function(x){which(x>0)})
  ####more processing to make nimble behave with n.cov=1 and/or n.rep=1
  
  G.obs=data$G.obs
  G.obs.NA.indicator=is.na(G.obs)
  G.obs[is.na(G.obs)]=9999999 #
  if(n.rep==1){
    warning("Padding G.obs (and NA indicator) along n.rep dimension to maintain array structure for Nimble. 2nd rep is all NA.")
    library(abind)
    G.obs=abind(G.obs,array(NA,dim=c(n.samples,n.cov,1)),along=3)
    G.obs.NA.indicator=abind(G.obs.NA.indicator,array(TRUE,dim=c(n.samples,n.cov,1)),along=3)
    n.rep=2 #need this if n.cov also 1
  }
  if(n.cov==1){
    warning("Padding G.obs (and NA indicator, and thetaArray) along n.cov dimension to maintain array structure for Nimble. 2nd cov is all NA.")
    warning("G.true and G.latent also padded. Ignore estimates for 2nd cov.")
    library(abind)
    G.obs=abind(G.obs,array(NA,dim=c(n.samples,1,n.rep)),along=2)
    G.obs.NA.indicator=abind(G.obs.NA.indicator,array(TRUE,dim=c(n.samples,1,n.rep)),along=2)
    G.true=cbind(G.true,rep(1,M))
    G.latent=cbind(G.latent,TRUE)
    thetaArray=abind(thetaArray,array(NA,dim=c(1,max(n.levels),max(n.levels))),along=1)
  }
  return(list(y.true=y.true2D,z=z,G.true=G.true,s=s,ID=ID,n.samples=n.samples,
              xlim=xlim,ylim=ylim,this.j=this.j,G.latent=G.latent,G.obs=G.obs,G.obs.NA.indicator=G.obs.NA.indicator,
              thetaArray=thetaArray))
}
