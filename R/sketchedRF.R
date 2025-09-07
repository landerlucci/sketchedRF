#' Sketched Random Forest
#'
#' Sketched Random Forest, proposed by Anderlucci, Montanari (2025), yields a forest of classification trees trained on multiple sets of data perturbed via matrix sketching.
#' @param x Data matrix of predictors.
#' @param y Vector of responses.
#' @param nk Vector of desired classes' size after sketching. By default, it is set as equal for all the classes, and equal to size of the largest class in the training set.
#' @param sketch Type of sketching matrix, one among \code{"Gaussian"},\code{"CW"}, \code{"Hadamard"}.
#' @param ntree no. of trees in the forest. Default is 500.
#' @param OOB Return sketched out-of-bag error estimate. Default is TRUE.
#' @param nk.OOB Desired classes' size for the out-of-bag sample. By default it is equal to the training set class distribution.
#' @param seed Set the seed to obtain reproducible results.
#' @return An object of class \code{sketchedRF}, which is a list with the following components:
  #' \item{call}{the original call to randomForest}
  #' \item{ntree}{number of trees grown}
  #' \item{mtry}{number of predictors sampled for spliting at each node.}
  #' \item{sketch}{type ok employed sketching matrices}
  #' \item{sOOB.err}{the sketched out-of-bag estimate of  error (when \code{OOB=TRUE}).}
#' @references L. Anderlucci, A. Montanari, M. Ferracin (2025). \emph{Perturbing data to address dataset shift in metastatic cancer classification}.
#' @seealso \link{predict.sketchedRF}, \link{varImp}, \link{varImpPlot.sketchedRF}
#' @examples
#' train<-c(1:45,51:95,101:145)
#' test<-c(46:50,96:100,146:150)
#' out.sRF<-sketchedRF(x=iris[train,-5],y=iris[train,5],nk=NULL,sketch="Gaussian",ntree=500,OOB=TRUE)
# 'out.sRF
#' @importFrom randomForest randomForest
#' @importFrom randomForest importance
#' @export

sketchedRF<-function(x=x,y=y,nk=NULL,sketch="Gaussian",ntree=500,OOB=TRUE,nk.OOB=NULL,seed=1){
  df.orig<-data.frame(y=as.factor(y),x)
  no.cl<-length(unique(y)) #no. classes
  if (is.null(nk)) {
    n0<-max(table(df.orig$y)) # size of the largest class
    nk<-rep(n0,no.cl) # final class sizes
  }
  if (!(sketch %in% c("Gaussian","gaussian","cw","CW","hada","Hada","hadamard","Hadamard"))) stop("If provided, argument 'scketch' must be one among 'Gaussian','CW' and 'Hadamard'")
  else {
    rf1.out<-list(Individual=list(),Sketch=sketch,x=x,y=y)

    if (OOB) {
      if (is.null(nk.OOB)) nk.OOB<-table(y)
      sk.OOB<-mMaSk(x,y,K=nk.OOB,type=sketch)
      votes.OOB<-matrix(NA,nrow(sk.OOB$x),ntree)
    }

    set.seed(seed)
    for(j in 1:ntree){
      sk.data<-mMaSk(x,y,K=nk,type=sketch)
      df.data<-data.frame(y=as.factor(sk.data$y),sk.data$x)
      rf1.out$Individual[[j]]<-randomForest(y~.,data=df.data,subset=1:nrow(df.data),replace=F,
                                 sampsize=nrow(df.data),ntree=1,keep.forest=T,importance=T)
      if (OOB) votes.OOB[,j]<-predict(rf1.out$Individual[[j]], newdata=sk.OOB$x,type="response")
    }
    if (OOB) {
      y.OOB<-apply(votes.OOB,1,getmode)
      ystar.OOB<-as.factor(sk.OOB$y)
      yhatOOB<-levels(ystar.OOB)[y.OOB]
      rf1.out$sOOB.err<-round(mean(ystar.OOB!=yhatOOB)*100,2)
      cat('sketched Random Forest with ',ntree,' trees. \n')
      cat(paste0('OOB estimate of  error rate: ',rf1.out$sOOB.err,'%.\n'))
    }

    class(rf1.out)<-"sketchedRF"
    }
return(rf1.out)
  }



