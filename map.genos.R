map.genos=function(x,unique.genos){
  n.cov=length(unique.genos)
  if(is.matrix(x)){ #is x a matrix?
    if(dim(x)[2]!=n.cov)stop("if x is a matrix, loci should be along the columns")
    x2=x*NA
    n.rows=nrow(x)
    for(j in 1:n.cov){
      x2[,j]=apply(unique.genos[[j]][x[,j],],1,paste,collapse="-")
    }
  }else if(is.array(x)){
    stop("Function does not work for arrays, yet")
  }else{
    x2=x*NA
    for(j in 1:n.cov){
      x2[j]=paste(unique.genos[[j]][x[j],],collapse="-")
    }
  }
  return(x2)
}