#' Variable Importance Plot
#'
#' Dotchart of variable importance as measured by a sketched Random Forest.
#' @param imp_obj an object of class sketchedRF, as that created by the function \code{importance}.
#' @param sort Should the variables be sorted in decreasing order of importance?
#' @param n.var	How many variables to show? (Ignored if sort=FALSE.)
#' @param main -plot title.
#' @param ...	Other graphical parameters to be passed on to dotchart.
#' @return Invisibly, the importance of the variables that were plotted.
#' @seealso \link{sketchedRF}, \link{predict.sketchedRF}, \link{varImp}
#' @examples
#' train<-c(1:45,51:95,101:145)
#' test<-c(46:50,96:100,146:150)
#' out.sRF<-sketchedRF(x=iris[train,-5],y=iris[train,5],nk=NULL,sketch="Gaussian",ntree=500)
#' iris_import<-varImp(out.sRF,type="all")
#' varImpPlot.sketchedRF(iris_import)
#' @importFrom graphics par
#' @importFrom graphics dotchart
#' @export

varImpPlot.sketchedRF<-function (imp_obj, sort = TRUE, n.var = min(30, nrow(imp_obj)),
                      main = deparse(substitute(imp_obj)),...)
{
  if (!inherits(imp_obj, "sketchedRF"))
                        stop("This function only works for objects of class `sketchedRF'")
  op <- par(mfrow = c(1, 2), mar = c(4, 5, 4, 1), mgp = c(2,
                                                        0.8, 0), oma = c(0, 0, 2, 0), no.readonly = TRUE)
  on.exit(par(op))
  for (i in 1:ncol(imp_obj)) {
    ord <- rev(order(imp_obj[, i], decreasing = TRUE)[1:n.var])
    xmin <-  min(imp_obj[ord, i])
    dotchart(imp_obj[ord, i], xlab = colnames(imp_obj)[i], ylab = "",
           main = NULL, xlim = c(xmin, max(imp_obj[, i])))
    }
  invisible(imp_obj)
  }
