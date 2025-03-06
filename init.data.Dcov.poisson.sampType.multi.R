init.data.Dcov.poisson.sampType.multi <- function(data=NA,M=NA,inits=inits){
  library(abind)
  #initialize sessions one at a time
  N.session <- length(data)
  init.sessions <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.sessions[[g]] <- init.data.Dcov.poisson.sampType(data=data[[g]],M=M[g],inits=inits.use)
  }
  
  ##Put data and inits into nimble-friendly structures (mostly ragged arrays)
  ##Preallocate structures so we can loop over sessions and fill them in##
  
  #typical SCR structure
  N.session <- length(data)
  if(length(M)!=N.session)stop("M and data must be of length 'N.session'")
  M.max <- max(M)
  J <- unlist(lapply(data,function(x){nrow(x$X)}))
  J.max <- max(J)
  X <- array(NA,dim=c(N.session,J.max,2))
  K <- rep(NA,N.session)
  K1D <- matrix(0,N.session,J.max)
  
  #genotype structures
  #these are not session-specific, pulling from first session
  n.cov <- data[[1]]$n.cov
  n.levels <- data[[1]]$n.levels
  n.levels.max <- max(n.levels)
  IDcovs <- data[[1]]$IDlist$IDcovs
  samp.levels <- nrow(inits$p.geno.het)
  #these are session-specific
  n.samples <- unlist(lapply(data,function(x){length(x$this.j)}))
  n.samples.max <- max(n.samples)
  n.rep <- unlist(lapply(data,function(x){dim(x$G.obs)[3]}))
  n.rep.max <- max(n.rep)
  if(n.rep.max==1){ #need to pad this if no reps anywhere
    n.rep.max <- 2
  }
  idx <- which(n.rep==1)
  if(length(idx)>0){
    n.rep[idx] <- 2 #bc padded
  }
  
  this.j <- this.k <- samp.type <- matrix(NA,N.session,n.samples.max)
  G.obs <- array(NA,dim=c(N.session,n.samples.max,n.cov,n.rep.max))
  #init structures
  y.true <- array(NA,dim=c(N.session,M.max,J.max))
  z <- matrix(NA,N.session,M.max)
  G.true <- G.latent <- array(NA,dim=c(N.session,M.max,n.cov))
  s <- array(NA,dim=c(N.session,M.max,2))
  ID <- matrix(NA,N.session,n.samples.max)
  G.obs.NA.indicator <- array(NA,dim=c(N.session,n.samples.max,n.cov,n.rep.max))
  thetaArray <- array(NA,dim=c(N.session,n.cov,samp.levels,n.levels.max,n.levels.max))
  
  n.cells <- unlist(lapply(data,function(x){x$n.cells}))
  n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
  n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
  n.cells.max <- max(n.cells)
  n.cells.x.max <- max(n.cells.x)
  n.cells.y.max <- max(n.cells.y)
  res <- unlist(lapply(data,function(x){x$res}))
  cellArea <- res^2
  x.vals <- matrix(NA,N.session,n.cells.x.max)
  y.vals <- matrix(NA,N.session,n.cells.y.max)
  dSS <- array(NA,dim=c(N.session,n.cells.max,2))
  InSS <- array(0,dim=c(N.session,n.cells.max))
  D.cov <- array(NA,dim=c(N.session,n.cells.max))
  cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
  xlim <- ylim <- matrix(NA,N.session,2)
  
  #Now fill in for each session
  for(g in 1:N.session){
    X[g,1:J[g],1:2] <- data[[g]]$X
    K[g] <- data[[g]]$K
    # K1D[g,1:J[g]] <- rep(K[g],J[g])
    K1D[g,1:J[g]] <- data[[g]]$K1D
    #genotyping structures - data
    this.j[g,1:n.samples[g]] <- data[[g]]$this.j
    this.k[g,1:n.samples[g]] <- data[[g]]$this.k
    samp.type[g,1:n.samples[g]] <- data[[g]]$samp.type
    #genotyping structures - inits
    y.true[g,1:M[g],1:J[g]] <- init.sessions[[g]]$y.true
    z[g,1:M[g]] <- init.sessions[[g]]$z
    G.true[g,1:M[g],] <- init.sessions[[g]]$G.true
    G.latent[g,1:M[g],] <- init.sessions[[g]]$G.latent
    G.obs.NA.indicator[g,1:n.samples[g],,1:n.rep[g]] <- init.sessions[[g]]$G.obs.NA.indicator
    G.obs[g,1:n.samples[g],,1:n.rep[g]] <- init.sessions[[g]]$G.obs
    s[g,1:M[g],] <- init.sessions[[g]]$s
    ID[g,1:n.samples[g]] <- init.sessions[[g]]$ID
    thetaArray[g,,,,] <- init.sessions[[g]]$thetaArray
    #Dcov structures
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
    x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
    y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
    dSS[g,1:n.cells[g],] <- data[[g]]$dSS
    InSS[g,1:n.cells[g]] <- data[[g]]$InSS
    D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
    cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
  }
  dummy.data <- matrix(0,N.session,M.max) #dummy data not used, doesn't really matter what the values are
  
  return(list(y.true=y.true,z=z,G.true=G.true,s=s,ID=ID,n.samples=n.samples,K1D=K1D,X=X,samp.type=samp.type,
              xlim=xlim,ylim=ylim,this.j=this.j,G.latent=G.latent,G.obs=G.obs,G.obs.NA.indicator=G.obs.NA.indicator,
              thetaArray=thetaArray,res=res,cellArea=cellArea,x.vals=x.vals,
              y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data))

}
