#' divided difference
#'
#' This function is to calculate single value of regressor
#'
#' @param t input x value
#' @param x knots
#' @return regressor value
#' @export




spval <- divided.difference <- function(t, x) {
  # dividived difference with order k = length(x) - 1
  # t is scalar, x is vector of length k+1
  # x should be at least length 2, and in non decreasing order
  # allow duplicated values in x
  # return scalar

  k <- length(x) - 1
  end_knot <- max(x)
  phi <- sapply(x, positive.polynomial, end_knot = end_knot,t=t, k=k)

  for (i in 1:k) {
    # length(phi_new) is always length(phi)-1
    phi_new <- rep(NA, k+1-i)
    for (j in 1:(k+1-i)) {
      if (x[j]==x[j+i]) {
        phi_new[j] <- abs(positive.polynomial(t, x[j], k, end_knot ,i))/factorial(i)
      } else {
        phi_new[j] <- abs((phi[j+1] - phi[j])/(x[j+i] - x[j]))
      }
    }
    phi <- phi_new
  }
  return((end_knot-x[1])* phi)
}
