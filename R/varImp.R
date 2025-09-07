#' Extract variable importance measure
#'
#' This is the function to compute variable importance measures as produced by sketched random Forest.
#' @param object an object of class skecthedRF, as that created by the function sketchedRF.
#' @param type one of 1, 2, "all", indicating the variable importance measure: \code{1}=mean decrease in accuracy, \code{2}=mean decrease in node impurity, \code{"all"} both of them.
#' @param nk Vector of desired classes' sizes for the sketched OOB sample. By default is the same of the training set.
#' @param seed Set the random seed for the generation of sketched OOB data.
#' @details Here are the definitions of the variable importance measures. The first measure is computed from permuting OOB data: For each tree, the prediction error on the out-of-bag portion of the data is recorded (error rate for classification). Then the same is done after permuting each predictor variable. The difference between the two are then averaged over all trees, and normalized by the standard deviation of the differences. The second measure is the total decrease in node impurities from splitting on the variable, averaged over all trees. The node impurity is measured by the Gini index.
#' @return A matrix of importance measure, one row for each predictor variable. The column(s) are different importance measures.
#' @references L. Anderlucci, A. Montanari, M. Ferracin (2025). \emph{Perturbing data to address dataset shift in metastatic cancer classification}.
#' @seealso \link{sketchedRF}, \link{predict.sketchedRF}, \link{varImpPlot.sketchedRF}
#' @examples
#' train<-c(1:45,51:95,101:145)
#' test<-c(46:50,96:100,146:150)
#' out.sRF<-sketchedRF(x=iris[train,-5],y=iris[train,5],nk=NULL,sketch="Gaussian",ntree=500)
#' varImp(out.sRF,type="all")
#' @importFrom stats predict
#' @importFrom stats sd
#' @export
varImp<-function(object, type="all", nk=NULL,seed=1){
  ntree<-length(object$Individual)
  x<-object$x
  y<-object$y
  imp<-matrix(NA,ncol(x),ntree,dimnames=list(colnames(x),NULL))
  if (is.factor(y)) y<-as.numeric(y)
  for (h in 1:ntree){ # h is the generic tree
    imp[,h]<-object$Individual[[h]]$importance[,ncol(object$Individual[[h]]$importance)]
  }
  varImp.gini<-apply(imp,1,mean)

 if (!(type==2)){
    if (is.null(nk)) nk<-table(y)
    sketch<-object$Sketch

    # Creation of the sketched OOB sample
    set.seed(seed)
    sk.data<-mMaSk(x,y,K=nk,type=sketch)
    sk.data.p<-sk.data$x
    votes.OOB<-votes.OOBp<-matrix(NA,nrow(sk.data$x),ntree)
    err.OOB<-rep(NA,ntree)
    varImp.OOB<-rep(NA,ncol(x))
    err.OOBp<-matrix(NA,ncol(x),ntree)

    # Accuracy of the OOB sample
    for (h in 1:ntree){ # h is the generic tree
      votes.OOB[,h]<-predict(object$Individual[[h]], newdata=sk.data$x,type="response")
      err.OOB[h]<-mean(votes.OOB[,h]!=sk.data$y)
    }

    # Accuracy of the OOB sample when permutations are computed
    for (j in 1:ncol(x)){
      # Permutation of the values of variable j
      sk.data.p[,j]<-sample(sk.data$x[,j],nrow(sk.data$x),replace=F)

      for (h in 1:ntree){
        # prediction of the permuted OOB for tree h
        votes.OOBp[,h]<-predict(object$Individual[[h]], newdata=sk.data.p,type="response")
        err.OOBp[j,h]<-mean(votes.OOBp[,h]!=sk.data$y) # classification rate for tree h
      }

      # The decrease in accuracy is averaged over all trees
      varImp.OOB[j]<-mean(err.OOBp[j,]-err.OOB)/sd(err.OOBp[j,]-err.OOB)
    }
  }

  if (type==2) {
    varImp<-matrix(varImp.gini,ncol=1)
    colnames(varImp)<-c("MeanDecreaseGini")
  } else {
    if (type==1) {
      varImp<-matrix(varImp.OOB,ncol=1)
      colnames(varImp)<-c("MeanDecreaseAccuracy")
    }
     else {
       varImp<-cbind(varImp.OOB,varImp.gini)
       colnames(varImp)<-c("MeanDecreaseAccuracy","MeanDecreaseGini")
     }
  }

  rownames(varImp)<-colnames(x)
  class(varImp)<-"sketchedRF"
  return(varImp)
  }
