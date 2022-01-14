build.genos=function(unique.genos){
  if(!is.list(unique.genos))stop("unique.genos should be a list of length n.cov (number of genetic loci)")
  n.cov=length(unique.genos)
  n.levels=unlist(lapply(unique.genos,nrow))
  zygotype=list(n.cov)
  ptype=list(n.cov)
  for(m in 1:n.cov){
    zygotype[[m]]=rep(NA,n.levels[m])
    ptype[[m]]=matrix(NA,n.levels[m],n.levels[m])
    for(l in 1:n.levels[m]){
      if(unique.genos[[m]][l,1]==unique.genos[[m]][l,2]){
        zygotype[[m]][l]=1 #homozygote
      }else{
        zygotype[[m]][l]=2 #heterozygote
      }
      for(l2 in 1:n.levels[m]){
        matches=sum(unique.genos[[m]][l,]==unique.genos[[m]][l2,])
        if(matches==2){#correct classification
          ptype[[m]][l,l2]=1 
        }else if(matches==1&
                 length(unique(unique.genos[[m]][l,]))==2&
                 length(unique(unique.genos[[m]][l2,]))==1){#allelic dropout.
          ptype[[m]][l,l2]=2
        }else{#false allele
          ptype[[m]][l,l2]=3
        }
      }
    }
  }
  #make a version of ptype that is a ragged array for nimble
  max.levels=max(n.levels)
  ptypeArray=array(NA,dim=c(n.cov,max.levels,max.levels))
  for(m in 1:n.cov){
    ptypeArray[m,1:n.levels[m],1:n.levels[m]]=ptype[[m]]
  }
  return(list(zygotype=zygotype,ptype=ptype,ptypeArray=ptypeArray))
}