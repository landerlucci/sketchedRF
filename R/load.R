#' Functions to be loaded
#'
#' Internal function that run behind sketchedRF
#' @param x data matrix
#' @param y vector of class labels
#' @param K vector of desired classes' size after sketching.
#' @param type type of matrix sketching; one among "Gaussian" (default), "CW" (Clarkson-Woodruff), "Hadamard".
#' @noRd

# Returned output: a list
## x: the sketched data matrix;
## y: the vector of class labels for the sketched data;
## K: vector of classes' size after sketching.;
## type: type of matrix sketching.

mMaSk<-function(x,y,K=NULL,type="Gaussian"){
  classi<-names(table(y))
  nk<-table(y)
  n<-sum(nk)
  k<-length(classi)
  p<-ncol(x)
  xx<-zz<-mu<-S<-X.tilde.star<-list()

  if (is.null(K)) K<-nk
  for (i in 1:k){
    xx[[i]]<-x[y==classi[i],]
    zz[[i]]<-scale(xx[[i]],T,F)
    mu[[i]]<-colMeans(xx[[i]])

    if (type %in% c('Gaussian','gaussian', 'g')){
      S[[i]]<-matrix(stats::rnorm(prod(nk[i],K[i]),sd=1/sqrt(K[i])),K[i],nk[i])
    }  else if (type %in% c('Clarkson-Woodruff','CW','cw')){
      S[[i]]<-clark.sk(nk[i],K[i])
    } else if (type %in% c('Hadamard','Hada','hadamard','hada')){
      S[[i]]<-hada.sk(nk[i],K[i])
      zz[[i]]<-rbind(zz[[i]],matrix(0,nrow=zero.row(nk[i]),ncol=p))
    }

    X.tilde.star[[i]]=sqrt(K[i]/nk[i])*(S[[i]]%*%zz[[i]])+matrix(mu[[i]],K[i],ncol=p,byrow=T)
  }
  xx.tilde.star<-do.call("rbind", X.tilde.star)
  y.star<-rep(classi,K)

  return(list(x=xx.tilde.star,y=y.star,nk=K,type=type))
}

clark.sk<-function(n,k){
#  require(extraDistr)
  I<-matrix(0,k,n)
  indice<-sample(1:k,n,replace=ifelse(n<k,FALSE,TRUE))
  for (i in 1:n) I[indice[i],i]<-extraDistr::rsign(1)
  return(I)
}


hada.sk<-function(n,k){
#  require(extraDistr)
#  require(pracma)
  if (zero.row(n)!=0) n<-n+zero.row(n)
  D<-extraDistr::rsign(n)
  quali.neg<-which(D<0)
  quali<-sample(1:n,k,replace=ifelse(n<k,TRUE,FALSE))
  IHD<-pracma::hadamard(n)[quali,]
  IHD[,quali.neg]<--IHD[,quali.neg]
  return(1/sqrt(k)*IHD)
}

# Function that returns the no. of zero rows to add in order to complete the Hadamard matrix
zero.row<-function(n){
  if (n %in% c(1,3,5,6,9,10))
    return(12-n)
  esponente=ceiling(log2(c(n,n/12,n/20)))
  xx=c(1,12,20)*2^esponente
  return(min(xx-n))
}

# Multiclass OverSampling
mOS<-function(x,y,K=NULL){
  classi<-names(table(y))
  nk<-table(y)
  n<-sum(nk)
  k<-length(classi)
  p<-ncol(x)

  if (is.null(K)) K<-nk
  X.os<-list()
  for (i in 1:k){
    X.os[[i]]<-x[sample(which(y==classi[i]),K[i],replace = T),]
  }
  x.os<-do.call("rbind", X.os)
  y.star<-rep(classi,K)

  return(list(x=x.os,y=y.star,K=K))
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode3<-function(v){
  uniqv <- unique(v)
  ni<-tabulate(match(v, uniqv))
  uniqv[order(ni,decreasing = T)[1:3]]
}


