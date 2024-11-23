e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

get.area <- function (X, buff){
  N.session <- length(X)
  area <- rep(NA,N.session)
  for(g in 1:N.session){
    area[g] <- diff(range(X[[g]][,1])+c(-buff[g],buff[g]))*diff(range(X[[g]][,2])+c(-buff[g],buff[g]))
  }
  return(area)
}

sim.genoSPIM.sampType.multi <-
  function(N.session=NA,N=NA,lam0=NA,sigma=NA,theta.d=NA,K=NA,p.geno.het=NA,
           p.geno.hom=NA,X=NA,buff=NA,n.cov=NA,n.rep=NA,p.sampType=NA,
           pID=NA,gamma=NA,IDcovs=NA,ptype=NA,obstype="poisson"){
    data <- vector("list",N.session)
    if(obstype=="poisson"){
      theta <- rep(NA,N.session)
    }
    for(g in 1:N.session){
      data[[g]] <- sim.genoSPIM.sampType(N=N[g],lam0=lam0[g],sigma=sigma[g],theta=theta[g],
                                      K=K[g],X=X[[g]],buff=buff[g],obstype=obstype,
                                      n.cov=n.cov,pID=pID,n.rep=n.rep[g],
                                      p.geno.hom=p.geno.hom[[g]],p.geno.het=p.geno.het[[g]],
                                      p.sampType=p.sampType,
                                      gamma=gamma[[g]],IDcovs=IDcovs,ptype=ptype)
    }
    return(data)
  }

