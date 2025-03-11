GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- (s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts <- numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1),G.latent=double(2)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

#build theta. Not as efficient as possible because we rebuild both
#homozygotes and heterozygotes every time p.geno.het or p.geno.hom elements are updated.
getTheta <- nimbleFunction(
  run = function(ptype = double(3), p.geno.het = double(1), p.geno.hom = double(1), n.levels = double(1)) {
    returnType(double(3))
    n.cov <- nimDim(ptype)[1]
    max.levels <- nimDim(ptype)[2]
    thetaArray <- array(NA,dim=c(n.cov,max.levels,max.levels))
    for(m in 1:n.cov){
      for(l in 1:n.levels[m]){
        if(!any(ptype[m,l,1:n.levels[m]]==2)){#homozygote bc no possible allelic dropout
          tmp3 <- (1/sum(ptype[m,l,1:n.levels[m]]==3))
          for(l2 in 1:n.levels[m]){
            if(ptype[m,l,l2]==3){#false allele
              thetaArray[m,l,l2] <- p.geno.hom[2]*tmp3
            }else{#correct
              thetaArray[m,l,l2] <- p.geno.hom[1]
            }
          }
        }else{ #heterozygote
          tmp2 <- (1/sum(ptype[m,l,1:n.levels[m]]==2))
          tmp3 <- (1/sum(ptype[m,l,1:n.levels[m]]==3))
          for(l2 in 1:n.levels[m]){
            if(ptype[m,l,l2]==3){#false allele
              thetaArray[m,l,l2] <- p.geno.het[3]*tmp3
            }else if(ptype[m,l,l2]==2){#allelic dropout
              thetaArray[m,l,l2] <- p.geno.het[2]*tmp2
            }else{#correct
              thetaArray[m,l,l2] <- p.geno.het[1]
            }
          }
        }
      }
    }
    return(thetaArray)
  }
)

dcat2 <- nimbleFunction(
  run = function(x = double(1), theta = double(1), n.rep = integer(0),
                 na.ind=double(1), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    for(rep in 1:n.rep){
      if(na.ind[rep]==FALSE){
        logProb <- logProb + log(theta[x[rep]])
      }
    }
    return(logProb)
  }
)
#dummy rng to make nimble happy, not used
rcat2 <- nimbleFunction(
  run = function(n = integer(0), theta = double(1), n.rep = integer(0),
                 na.ind=double(1)) {
    returnType(double(1))
    out <- numeric(n.rep,value=999)
    return(out)
  }
)


#------------------------------------------------------------------
# Custom sampler to update G.true, using information in G.latent to determine proposal distribution
#------------------------------------------------------------------
#Metropolis-Hastings update here allows other parameters to be a function of G.true
#Unsure if this is the most efficient way to calculate the likelihood for the proposal.
#I assume nimble is attempting to update the likelihood of every G.obs[i,m] and skipping
#the samples not involved. But we know which ones are involved. "these.samps" below.
#Anyways, perhaps an inefficiency here.
GSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    n.levels <- control$n.levels
    n.rep <- control$n.rep
    na.ind <- control$na.ind
    calcNodes <- model$getDependencies(target)
    i <- control$i
    m <- control$m
    samp.type <- control$samp.type
  },
  run = function() {
    G.probs <- model$gammaMat[m,1:n.levels[m]] #pull out genotype frequencies
    if(model$G.latent[i,m]==FALSE){ #these individual-loci have samples allocated to them currently
      #build proposal distribution using genotype frequency likelihood and classification likelihood.
      these.samps <- which(model$ID==i)
      G.obs <- model$G.obs
      error.probs=rep(1,n.levels[m]) #error probs|category level
      for(i2 in 1:length(these.samps)){
        for(obs in 1:n.rep){
          if(na.ind[these.samps[i2],obs]==FALSE){ #if observed
            error.probs <- error.probs*model$theta[m,samp.type[these.samps[i2]],1:n.levels[m],G.obs[these.samps[i2],m,obs]]
          }
        }
      }
      G.probs <- G.probs*error.probs
      G.probs <- G.probs/sum(G.probs)
    }else{ #these individual-loci do not have samples allocated to them currently
      #build proposal distribution using only genotype frequency likelihood. This is the
      #full conditional if no other parameters depend on G.true
      G.probs <- G.probs/sum(G.probs)
    }
    #MH step
    model.lp.initial <- model$getLogProb(calcNodes) #initial logProb
    prop.back <- G.probs[model$G.true[i,m]] #backwards proposal prob
    G.prop <- rcat(1,G.probs[1:n.levels[m]])
    prop.for <- G.probs[G.prop] #forwards proposal prob
    model$G.true[i,m] <<- G.prop #store in model
    model.lp.proposed <- model$calculate(calcNodes)#proposed logProb
    log_MH_ratio <- (model.lp.proposed+log(prop.back)) - (model.lp.initial+log(prop.for))
    # log_MH_ratio
    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

#Alternative G.true sampler that uses much less RAM and is faster.
#However! No parameters other than G.obs can depend on G.true with this
#version unless you add their likelihoods to the MH update.
GSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    n.cov <- control$n.cov
    n.levels <- control$n.levels
    n.samples <- control$n.samples
    n.rep <- control$n.rep
    na.ind <- control$na.ind
    samp.type <- control$samp.type
    G.true.nodes <- control$G.true.nodes
    G.obs.nodes <- control$G.obs.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    node.idx <- 1 #this is used to keep up with the correct nodes in G.true.nodes which MUST BE 1D, for each i and m
    #must loop over loci first, then individuals for this to work correctly.
    for(m in 1:n.cov){
      G.obs.m.idx <- seq((m-1)*n.samples+1,(m-1)*n.samples+n.samples,1) #get index for cov m for all samples
      for(i in 1:M){
        G.probs <- model$gammaMat[m,1:n.levels[m]] #pull out genotype frequencies
        if(model$G.latent[i,m]==FALSE){ #these individual-loci have samples allocated to them currently
          #must use MH
          #build proposal distribution using genotype frequency likelihood and classification likelihood.
          these.samps <- which(model$ID==i)
          G.obs.idx <- G.obs.m.idx[these.samps] #pull out nodes for these samples at cov m
          n.these.samps <- length(these.samps)
          G.obs <- model$G.obs
          error.probs <- rep(1,n.levels[m]) #error probs|category level
          for(i2 in 1:length(these.samps)){
            for(obs in 1:n.rep){
              if(na.ind[these.samps[i2],m,obs]==FALSE){ #if observed
                error.probs <- error.probs*model$theta[m,samp.type[these.samps[i2]],1:n.levels[m],G.obs[these.samps[i2],m,obs]]
              }
            }
          }
          G.probs <- G.probs*error.probs
          G.probs <- G.probs/sum(G.probs)
          #MH step
          G.true.lp.initial <- model$getLogProb(G.true.nodes[node.idx]) #initial logProb for G.true
          G.obs.lp.initial <- model$getLogProb(G.obs.nodes[G.obs.idx]) #initial logProb for G.obs
          prop.back <- G.probs[model$G.true[i,m]] #backwards proposal prob
          G.prop <- rcat(1,G.probs[1:n.levels[m]]) #proposal
          if(G.prop!=model$G.true[i,m]){ #we can skip this if we propose the current value
            prop.for <- G.probs[G.prop] #forwards proposal prob
            model$G.true[i,m] <<- G.prop #store in model
            G.true.lp.proposed <- model$calculate(G.true.nodes[node.idx])#proposed logProb for G.true
            G.obs.lp.proposed <- model$calculate(G.obs.nodes[G.obs.idx])#proposed logProb for G.true
            log_MH_ratio <- (G.true.lp.proposed+G.obs.lp.proposed+log(prop.back)) -
              (G.true.lp.initial+G.obs.lp.initial+log(prop.for))
            # log_MH_ratio
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["G.true",1][i,m] <<- model[["G.true"]][i,m]
            } else {
              model[["G.true"]][i,m] <<- mvSaved["G.true",1][i,m] #set back to init
              model$calculate(G.true.nodes[node.idx]) #set log prob back to init
              model$calculate(G.obs.nodes[G.obs.idx]) #set log prob back to init
            }
          }
        }else{ #these individual-loci do not have samples allocated to them currently
          #build proposal distribution using only genotype frequency likelihood. This is the
          #full conditional if no other parameters depend on G.true. So always accept
          G.probs=G.probs/sum(G.probs)
          G.prop <- rcat(1,G.probs[1:n.levels[m]])
          model$G.true[i,m] <<- G.prop #store in model
          model$calculate(G.true.nodes[node.idx]) #update G.true logprob. No G.obs logprob.
          mvSaved["G.true",1][i,m] <<- model[["G.true"]][i,m]
        }
        node.idx=node.idx+1 #increment
      }
    }
    #copy back to mySaved to update logProbs. should be done above, already though.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#------------------------------------------------------------------
# Customer sampler to update latent IDs, and associated arrays
#------------------------------------------------------------------
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    J <- control$J
    K1D <- control$K1D
    G.obs <- control$G.obs
    n.cov <- control$n.cov
    n.rep <- control$n.rep
    n.levels <- control$n.levels
    n.samples <- control$n.samples
    this.j <- control$this.j
    na.ind <- control$na.ind
    samp.type <- control$samp.type
    calcNodes <- model$getDependencies(c("y.true","G.obs"))
  },
  run = function() {
    G.true <- model$G.true
    y.true <- model$y.true
    ID.curr <- model$ID
    z <- model$z
    theta <- model$theta
    lam <- model$lam
    #Precalculate likelihoods
    ll.y <- matrix(0,nrow=M,ncol=J) #M x J most efficient here
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,1:J] <- dpois(y.true[i,1:J],K1D[1:J]*lam[i,1:J],log=TRUE)
      }
    }
    ll.y.cand <- ll.y
    ll.theta <- array(0,dim=c(n.samples,n.cov,n.rep))
    for(l in 1:n.samples){
      for(m in 1:n.cov){
        for(rep in 1:n.rep){
          if(na.ind[l,m,rep]==FALSE){
            ll.theta[l,m,rep] <- dcat(G.obs[l,m,rep],theta[m,samp.type[l],G.true[ID.curr[l],m],1:n.levels[m]],log=TRUE)
          }
        }
      }
    }
    ll.theta.cand <- ll.theta
    y.cand <- y.true
    for(l in 1:n.samples){
      #proposal distribution is combination of distance-based and genotype-based distributions
      dist.probs <- lam[,this.j[l]]*z
      #G.probs proportional to genotyping error likelihood over all loci and reps
      G.probs <- z #setting to z initializes to 1 if z=1, 0 otherwise
      for(i in 1:M){
        # if(z[i]==1){
        if(z[i]==1&dist.probs[i]>1e-10){ #skip if z=0 or s not close to focal s
          for(m in 1:n.cov){
            for(rep in 1:n.rep){
              if(na.ind[l,m,rep]==FALSE){ #if observed
                G.probs[i] <- G.probs[i]*theta[m,samp.type[l],G.true[i,m],G.obs[l,m,rep]]
              }
            }
          }
        }
      }
      total.probs <- dist.probs*G.probs
      total.probs <- total.probs/sum(total.probs) 
      ID.cand <- ID.curr
      ID.cand[l] <- rcat(1,prob=total.probs)
      if(ID.curr[l]!=ID.cand[l]){ #skip if propose same ID
        swapped <- c(ID.curr[l],ID.cand[l])#order swap.out then swap.in
        propprob <- total.probs[swapped[2]]
        backprob <- total.probs[swapped[1]]
        
        #update y.true
        y.cand[ID.curr[l],this.j[l]] <- y.true[ID.curr[l],this.j[l]]-1
        y.cand[ID.cand[l],this.j[l]] <- y.true[ID.cand[l],this.j[l]]+1
        #select sample to move proposal probabilities
        focalprob <- y.true[ID.curr[l],this.j[l]]/n.samples
        focalbackprob <- y.cand[ID.cand[l],this.j[l]]/n.samples
        
        ##update ll.y
        ll.y.cand[swapped,this.j[l]] <- dpois(y.cand[swapped,this.j[l]],K1D[this.j[l]]*model$lam[swapped,this.j[l]],log=TRUE)
        
        #update ll.theta
        for(m in 1:n.cov){
          for(rep in 1:n.rep){
            if(na.ind[l,m,rep]==FALSE){
              ll.theta.cand[l,m,rep] <- dcat(G.obs[l,m,rep],theta[m,samp.type[l],G.true[ID.cand[l],m],1:n.levels[m]],log=TRUE)
            }else{
              ll.theta.cand[l,m,rep] <- 0
            }
          }
        }
        if(runif(1)<exp((sum(ll.y.cand[swapped,this.j[l]])+sum(ll.theta.cand[l,,]))-
                        (sum(ll.y[swapped,this.j[l]])+sum(ll.theta[l,,])))*
           (backprob/propprob)*(focalbackprob/focalprob)){
          y.true[swapped,this.j[l]] <- y.cand[swapped,this.j[l]]
          ll.y[swapped,this.j[l]] <- ll.y.cand[swapped,this.j[l]]
          ll.theta[l,,] <- ll.theta.cand[l,,]
          ID.curr[l] <- ID.cand[l]
        }else{
          #set these back.
          y.cand[swapped,this.j[l]] <- y.true[swapped,this.j[l]]
        }
      }
    }
    
    #update G.latent after ID changes
    G.latent <- matrix(TRUE,nrow=M,ncol=n.cov)
    for(l in 1:n.samples){
      for(m in 1:n.cov){
        for(rep in 1:n.rep){
          if(na.ind[l,m,rep]==FALSE){ #if this sample is not latent
            G.latent[ID.curr[l],m] <- FALSE #then its ID is not latent
          }
        }
      }
    }
    #put everything back into the model$stuff after updating y.sight.true, y.sight.true.event
    model$y.true <<- y.true
    model$ID <<- ID.curr
    model$G.latent <<- G.latent
    # model.lp.proposed <-
    model$calculate(calcNodes) #update dependencies, likelihoods
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)