#' augknt
#'
#' This function is used to find knots
#'
#' @param boundary boundary of knots
#' @param nBspline number of b splines
#' @param order order of b spline function, indicates the degree n-1 of piecewise polynomial function
#' @export
#' @examples
#' boundary = c(0.1, 0.5)
#' nBspline = 10
#' order = 3
#' print(augknt(boundary, nBspline, order))
#'


augknt <- function(boundary, nBspline, order) {
  # return knots of length nBspline + order
  c(rep(boundary[1], order - 1),
    seq(boundary[1], boundary[2], length.out = nBspline-order+2),
    rep(boundary[2], order - 1))
}





