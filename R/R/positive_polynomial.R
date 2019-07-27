#' positive.polynomial
#'
#' This function works as truncated polynomial function
#'
#' @param t the input value
#' @param x knots value
#' @param k order
#' @param end_knot the largest value for knot
#' @param d the dth derivative if wanted
#' @return a number of positive polynomial result
#' @export
#' @examples
#' positive.polynomial(0.2,x = 0,k = 3,end_knot = 2,d=0)

positive.polynomial <- function (t, x, k,end_knot, d=0) {
  # positive polynomial of order k(degree k-1) cutting at x, building block of bspline
  # d is the dth derivative if wanted
  # t, x be scalar
  if (t <= x || d > k-1 ||t > end_knot) {
    return(0)
  } else if (d > 0) {
    return(exp(lfactorial(k-1)-lfactorial(k-d-1)) * (t-x)^(k-1-d))
  } else {
    return((t-x)^(k-1))
  }
}
