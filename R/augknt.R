#' augknt
#'
#' This function is to find knots
#'
#' @param boundary boundary of knots
#' @param nBspline number of b splines
#' @param order order of polynomial
#' @export


augknt <- function(boundary, nBspline, order) {
  # return knots of length nBspline + order
  c(rep(boundary[1], order - 1),
    seq(boundary[1], boundary[2], length.out = nBspline-order+2),
    rep(boundary[2], order - 1))
}





