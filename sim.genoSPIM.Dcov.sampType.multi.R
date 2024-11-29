sim.genoSPIM.Dcov.sampType.multi <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,sigma=NA,theta.d=NA,K=NA,p.geno.het=NA,
           p.geno.hom=NA,X=NA,xlim=NA,ylim=NA,n.cov=NA,n.rep=NA,p.sampType=NA,
           pID=NA,gamma=NA,IDcovs=NA,ptype=NA,obstype="poisson"){
    data <- vector("list",N.session)
    if(obstype=="poisson"){
      theta <- rep(NA,N.session)
    }
    for(g in 1:N.session){
      data[[g]] <- sim.genoSPIM.Dcov.sampType(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],
                                               InSS=InSS[[g]],xlim=xlim[g,],ylim=ylim[g,],res=res[g],
                                               lam0=lam0[g],sigma=sigma[g],theta=theta[g],
                                               K=K[g],X=X[[g]],obstype=obstype,
                                               n.cov=n.cov,pID=pID,n.rep=n.rep[g],
                                               p.geno.hom=p.geno.hom[[g]],p.geno.het=p.geno.het[[g]],
                                               p.sampType=p.sampType,gamma=gamma[[g]],IDcovs=IDcovs,ptype=ptype)
    }
    return(data)
  }

