#' get Preprocess
#'
#' This function preprocess the original data to get splines configuration and input data for getRegressor function.
#'
#' @param X The origianl input data.
#' @param L The maximum lag.
#' @param order order of b splines.
#' @param nBspline number of B splines.
#' @param ... defualt parameters,including
#'  \itemize{
#'   \item scaleFlag: boolean type,if true, standardize the columns of X, default to TRUE
#'   \item augkntFlag:if false, use linear regressor, defaul to TRUE
#' }
#' @import stats
#' @import methods
#' @return a list including
#'  \itemize{
#'   \item x: reformed data
#'   \item scale: standard deviation for original input
#'   \item feasibleBox: the region of X where the spline fit is very accurate
#'   \item knotBox: determine the knot range IN RESCALED SCALE, which is slightly larger than feasibleBox in quantile
#'   \item knots: knots
#'   \item spconfig: splines configuration for b splines
#' }
#' @export
#' @examples
#' utils.lag <- function(ts, lag = 1, pad = NA) {
#' # return the lagged version of a time series vector
#' return(c(rep(pad, lag), ts[1:(length(ts)-lag)]))
#' }

#'# pack configurations so that it is easier to bundle things together
#'
#' # Suppose that there is no historical data. In this situation, onlineModel has the same function as getSequentialNonlinearModel
#'Ex1 <- (function(N=2000, D=2, L=8, ifPlot = TRUE){
#'  err <- matrix(stats::rnorm(D*(N+L)), N+L, D)
#'  X <- err[ ,1]
#'  X <- cbind(X, 0.5 * utils.lag(X,1)^2 - 0.8 * utils.lag(X,7) + 0.2 * err[,2])
#'  if (ifPlot) {
#'    graphics::plot(X[,1], X[,2], pch = 16, cex = 0.5, col = "red")
#'  }
#'  list(N=N, D=D, L=L, X=X[-(1:L),], y=X[-(1:L),2])
#'})()
#'par <- c(order = 3, nBspline = 10)
#'result <- do.call(getPreprocess, c(Ex1,par))



getPreprocess <- function(X, L, order, nBspline, ...) {

  ###
  # X is original design matrix of dimension N by D, L is maximum lag
  # output design matrix x of dimension N by (D*L) by adding lag of length L for each column
  # output is a list of following:
  #   1) x of dim N by (D*L)
  #   2) scale, length D
  #   3) feasibleBox, used for plots
  #   4) knots, dim (nBspline+order) by D
  ###

  ########################
  # set default parameters
  dots <- list(...)

  # if true, standardize the columns of X, default to TRUE
  scaleFlag <- ifelse(methods::hasArg(scaleFlag), dots$scaleFlag, TRUE)
  # if false, use linear regressor
  augkntFlag <- ifelse(methods::hasArg(augkntFlag), dots$augkntFlag, TRUE)

  ########################
  N <- nrow(X)
  D <- ncol(X)

  # rescale (standardize) x
  x <- matrix(NA, nrow = D * L, ncol = N)

  if (scaleFlag) {
    scale <- apply(X, 2, stats::sd, na.rm = TRUE)
  } else {
    scale <- rep(1, D)
  }

  for (lag in 1:L) {
    # arrange the data matrix in row order: X(1,n-1)...X(D,n-1), ..., X(1,n-L)...X(D,n-L)
    x[((lag-1)*D+1):(lag*D), (L+1):N] <- t(X[(L+1-lag):(N-lag), ])/scale
  }
  # transpose x back to design matrix
  x <- t(x)


  ########################
  # determine the feasible box IN ORIGINAL SCALE---the region of X where the spline fit is very accurate
  feasibleBox <- apply(X, 2, stats::quantile, c(0.1,0.9), na.rm = TRUE)
  # determine the knot range IN RESCALED SCALE, which is slightly larger than feasibleBox in quantile
  knotBox = apply(x[,1:D], 2, stats::quantile, c(0.02,0.98), na.rm = TRUE)

  ########################
  # determine knots and spline basis
  if (augkntFlag) {
    knots <- apply(knotBox, 2, augknt, nBspline, order)
  } else {
    # linear predictor
    knots <- knotBox[c(1,2,2), ]
  }

  # determine the basis
  # there is no counterpart of spmak in R
  # we just store the spline configuration and generate it on the fly (slow though)
  # basis{j, (ell-1)*D + d} = spmak(knots(d, j:j+order), 1)
  spconfig = list(order = order, nBspline = nBspline, knots = knots,knotBox = knotBox)

  list(x=x, scale=scale, feasibleBox=feasibleBox,knotBox = knotBox, knots=knots, spconfig = spconfig)
}
