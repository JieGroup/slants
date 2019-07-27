#' divided difference
#'
#' This function is to calculate single value of regressor
#'
#' @param x input x value
#' @param knots knots
#' @return regressor value
#' @export
#' @examples
#' knots = c(0,0,0,0.1)
#' input = 0.02
#' result = spval(input,knots)
#'




spval <- divided.difference <- function(x, knots) {
  # dividived difference with order k = length(knots) - 1
  # x is scalar, knots is vector of length k+1
  # knots should be at least length 2, and in non decreasing order
  # allow duplicated values in knots
  # return scalar

  k <- length(knots) - 1
  end_knot <- max(knots)
  phi <- sapply(knots, positive.polynomial, end_knot = end_knot,x=x, k=k)

  for (i in 1:k) {
    # length(phi_new) is always length(phi)-1
    phi_new <- rep(NA, k+1-i)
    for (j in 1:(k+1-i)) {
      if (knots[j]==knots[j+i]) {
        phi_new[j] <- abs(positive.polynomial(x, knots[j], k, end_knot ,i))/factorial(i)
      } else {
        phi_new[j] <- abs((phi[j+1] - phi[j])/(knots[j+i] - knots[j]))
      }
    }
    phi <- phi_new
  }
  return((end_knot-knots[1])* phi)
}
