#' plotcoeff
#'
#' This function is to plot the coefficient
#'
#' @param beta optimal coefficients
#' @param knots knots that returned from getPreprocess function
#' @param nBspline the numbers of B splines that we need to fit model, also returned from getPreprocess function
#' @import graphics
#' @export
#' @examples
#' positive.polynomial(0.2,x = 0,k = 3,end_knot = 2,d=0)

plotcoeff <- function(beta,knots,nBspline){
  knots <- knots[,1]
  D <- length(knots)
  DL <- length(beta)/nBspline
  order <- 3
  x <- seq(-1,1,by=0.01)
  tmp <- NULL

  for (i in 1:(length(knots)-order)) {
    tmp <- cbind(tmp, sapply(seq(-1,1,by=0.01), spval, x=knots[i:(i+order)]))
  }
  result = matrix(nrow = DL,ncol = length(x))
  rownames(result) <- paste('lag',1:DL,sep='')
  for (i in 1:DL){
    betanew = beta[(nBspline *(i-1)+1):(nBspline*i)]
    result[i,] = tmp %*% betanew

  }

  graphics::plot(x, result[1,],ylim=c(-1,1),col = 1)
  for (i in 2:DL){
    graphics::lines(x, result[i,],col = i)
  }
  graphics::grid()
  graphics::legend("topright",legend=rownames(result),lty=c(1,2,3),col=1:DL,cex = 0.5)

}
