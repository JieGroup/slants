#' plotcoeff
#'
#' This function is to plot the coefficient
#'
#' @param beta optimal coefficients
#' @param knots first column of knots that returned from getPreprocess function, generated from original scale
#' @param nBspline the numbers of B splines that we need to fit model, also returned from getPreprocess function
#' @param order the order of b splines
#' @import graphics
#' @export
#' @examples
#' beta = c(rep(1,10),rep(0,140),rep(-1,10))
#' knots = c(0.1,0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.8,0.8,0.8)
#' nBspline = 10
#' order = 3
#' plotcoeff(beta,knots,nBspline,order)

plotcoeff <- function(beta,knots,nBspline,order){
  D <- length(knots)
  DL <- length(beta)/nBspline
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

  graphics::plot(x, result[1,],ylim=c(-1,1),col = 1,ylab = "Y")
  for (i in 2:DL){
    graphics::lines(x, result[i,],col = i)
  }
  graphics::grid()
  graphics::legend("topright",legend=rownames(result),lty=c(1,2,3),col=1:DL,cex = 0.5)

}
