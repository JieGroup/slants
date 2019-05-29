#' glasso_EM
#'
#' This function is to implement the em algorithm for group lasso
#'
#' @param beta_init initial value
#' @param A sufficient statistics, matrix
#' @param B sufficient statistics
#' @param ngroup number of groups
#' @param lambda the LASSO penalty
#' @param alpha2 the EM decomposition parameter, refer to the paper
#' @param ... default parameters
#' @import methods
#' @import graphics
#' @return a list of whether success and the updated coefficients
#' @export
#' @examples
#' beta = matrix(c(1,2,3))
#' B = matrix(c(1,2,3))
#' A = matrix(rep(c(1,2,3),3),nrow = 3,byrow = TRUE)
#' ngroup = 3
#' lambda = 0.05
#' alpha2 = 0.05
#' out = glasso_EM(beta,A,B,ngroup,lambda,alpha2)
#'





glasso_EM <- function(beta_init, A, B, ngroup, lambda, alpha2, ...) {
  ###
  # beta_init is the initial value
  # sufficient statistics A (matrix) and B (vector)
  # lambda is the LASSO penalty
  # ngroup is the number of groups, should be D*L in our setting
  # TODO: we currently only consider equal group size
  # alpha2 is the EM decomposition parameter, refer to the paper
  # return a list, consisting of
  #   1) whether success
  #   2) the updated coefficients
  ###

  ########################
  # set default parameters

  dots <- list(...)

  # K is the number of iterations, exponential convergence, small value is sufficient
  K <- ifelse(methods::hasArg(K), dots$K, 20)

  # tolerance.EM is the tolerance of EM convergence (whether explodes)
  tolerance.EM <- ifelse(methods::hasArg(tolerance.EM), dots$tolerance.EM, 100)

  # eps_EM is the convergence criteria of EM
  # first component is the absolute change threshold
  # second component is the relative change threshold
  eps_EM_abs <- ifelse(methods::hasArg(eps_EM_abs), dots$eps_EM_abs, 1e-3)
  eps_EM_rel <- ifelse(methods::hasArg(eps_EM_rel), dots$eps_EM_rel, 1e-2)

  # default: no print or plot
  ifPrint <- ifelse(methods::hasArg(ifPrint), dots$ifPrint, 0)

  ########################
  # total number of parameters and group size
  p <- length(beta_init)
  gsize <- p/ngroup

  # initialize the history container
  beta_history <- beta <- beta_init

  ########################
  # iterate between old beta -> compute r -> shrink r (new beta)
  for (k in 1:K) {
    beta_old <- beta
    r <- beta - alpha2 * A %*% beta + alpha2 * B
    for (j in 1:ngroup) {
      group <- ((j-1)*gsize+1):(j*gsize)
      temp <- alpha2 * lambda / sqrt(sum(r[group]^2))
      if (temp<1) {
        beta[group] <- (1-temp)*r[group]
      } else {
        beta[group] <- 0
      }
    }
    beta_history <- rbind(beta_history, beta)
    # check convergence
    chg <- sum(abs(beta - beta_old))
    if (chg < eps_EM_abs) {
      if (ifPrint) {
        print(sprintf("EM converges in %d iterations (abs)", k))
      }
      break
    } else if (sum(abs(beta_old))>0 && chg/sum(abs(beta_old)) < eps_EM_rel) {
      if (ifPrint) {
        print(sprintf("EM converges in %d iterations (rel)", k))
      }
      break
    }
  }

  ########################
  # post processing, check convergence
  if (sum(abs(beta_init))>0 && sum(abs(beta-beta_init)) > tolerance.EM * sum(abs(beta_init))) {
    success <- FALSE
    print(sprintf("EM blow up with rate %.4ef", sum(abs(beta-beta_init))/sum(abs(beta_init))))
  } else {
    success <- TRUE
  }

  if (ifPrint>1) {
    # plot the history of first component
    graphics::plot(beta_history[,1], main = "EM history of beta[1]", xlab = "Iteration", ylab = "beta[1]")
  }

  list(success = success, beta = beta)
}
