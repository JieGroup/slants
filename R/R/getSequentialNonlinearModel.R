#' getSequentialNonlinearModel
#'
#' This function uses splines to approximate a wide range of nonlinear functions and adaptive filtering to accommodate time varying data generating processes.
#'
#' @param y vector
#' @param x the reuslt from getPreprocess \code{x}
#' @param D columns of original X
#' @param L maximize lag number
#' @param forget_factor Comprises ideally decreasing weights in range [0,1] indicating the importance of past values (generally, less significance for older values)
#' @param lambda_init the initial value for lambda, LASSO penalty
#' @param tau2_init the initial value of EM decomposition parameter, refer to the paper
#' @param spconfig spline configuration
#' @param delta_multiplier multipler to move among channels
#' @param sparse_lambda_tolerance favor smaller lambda within spaTol_c tolerance
#' @param shrinkStepSize adjust lambda before they get too crazy, recommended to be equal to 'forget_factor' vector values in case of confusion
#' @param ... default parameters,including
#'  \itemize{
#'   \item testSize: the number of obs to average prediction error, default is 50
#'   \item tol_idle_tau2:increase tau2 if it idles for a while
#'   \item safeShrink_lambda : adjust lambda before they get too crazy
#'   \item tol_all0times: if coefficents are all 0 for too long time, shrink lambda
#'   \item shrinktau2: how to shrink tau2
#' }
#' @return a list including
#'  \itemize{
#'   \item preErr: historical predicted error, firt column indicates the first channel, etc.
#'   \item predict1 : historical predicted values, firt column indicates the first channel, etc.
#'   \item lambda_history : historical lambda channels, first channel is the smallest one and the third channel is the largest one.
#'   \item lambda_opt: optimized lambda, the second channel in lambda_history
#'   \item tau2_opt: optimized tau2
#'   \item beta_opt: best coefficient, the last row and the second channel related in beta_history
#'   \item beta_history: all stored coefficient.
#' }
#' @import stats
#' @import methods
#' @import graphics
#' @export
#' @description This model is used to fit sequential non linear time series data. Detailed information could be find in the paper \url{http://jding.org/jie-uploads/2018/11/slant.pdf}
#' @seealso \code{\link{glasso_EM}} for implementation of EM algorithm; \code{\link{getPreprocess}};\code{\link{getRegressor}}
#' @keywords  sequential non linear time series model
#' @examples
#' utils.lag <- function(ts, lag = 1, pad = NA) {
#' # return the lagged version of a time series vector
#' return(c(rep(pad, lag), ts[1:(length(ts)-lag)]))
#' }

#'# pack configurations so that it is easier to bundle things together

#'Ex1 <- (function(N=2000, D=2, L=8, ifPlot = TRUE){
#'  err <- matrix(stats::rnorm(D*(N+L)), N+L, D)
#'  X <- err[ ,1]
#'  X <- cbind(X, 0.5 * utils.lag(X,1)^2 - 0.8 * utils.lag(X,7) + 0.2 * err[,2])
#'  if (ifPlot) {
#'    graphics::plot(X[,1], X[,2], pch = 16, cex = 0.5, col = "red")
#'  }
#'  list(N=N, D=D, L=L, X=X[-(1:L),], y=X[-(1:L),2])
#'})()
#'
#'Ex1_algo <- (function(experiment_config){
#' ec <- experiment_config # short name
#' forget_factor <- 1/c(1:ec$N) # same as batch
#' shrinkStepSize <- 1/c(1:ec$N)
#' list(
#'   order = 3,
#'   nBspline = 10,
#'   forget_factor = forget_factor,
#'   shrinkStepSize = shrinkStepSize,
#'   # if the performance does not exceed this ratio, the sparser (larger lambda) is preferred
#'   sparse_lambda_tolerance = 1.01,
#'   delta_multiplier = 10^(0.4), # multipler to move among channels
#'  safeShrink_lambda = 10^(0.1), # to adjust lambda before they get too crazy
#'   lambda_init = 0.01,
#'   tau2_init = 0.05
#' )
#' })(Ex1)
#'
#'Ex1_algo <- c(Ex1_algo, do.call(getPreprocess, c(Ex1,Ex1_algo)))
#'par <- c(list(ifPrint=0,testSize=50),Ex1,Ex1_algo)
#'Ex1_result <- do.call(getSequentialNonlinearModel, par)
#'#========diagonistic========
#'graphics::plot(Ex1_result$lambda_opt,type = "l")
#'graphics::plot(Ex1_result$tau2_opt,type = "l",ylab = "Tao2")
#'graphics::plot(Ex1_result$preErr[,2],type = "l")
#'#get the historical opt beta, always the middle channel
#'plotcoeff(Ex1_result$beta_opt,Ex1_algo$knots[,1],Ex1_algo$spconfig$nBspline,Ex1_algo$spconfig$order)



getSequentialNonlinearModel <- function(y, x, D, L, forget_factor, lambda_init, tau2_init, spconfig, sparse_lambda_tolerance, shrinkStepSize, delta_multiplier, ...) {



  ########################
  # set default parameters


  dots <- list(...)

  sparse_lambda_tolerance <- ifelse(methods::hasArg(sparse_lambda_tolerance), sparse_lambda_tolerance, 1.01)

  lambda_init <- ifelse(methods::hasArg(lambda_init), lambda_init, 0.03)

  tau2_init <- ifelse(methods::hasArg(tau2_init), tau2_init, 0.05)

  # testSize is the number of obs to average prediction error, default is 50
  testSize <- ifelse(methods::hasArg(testSize), dots$testSize, 50)

  # increase tau2 if it idles for a while
  tol_idle_tau2 <- ifelse(methods::hasArg(tol_idle_tau2), dots$tol_idle_tau2, 50)

  # to adjust lambda before they get too crazy
  safeShrink_lambda <- ifelse(methods::hasArg(safeShrink_lambda), dots$safeShrink_lambda, 10^(0.4))

  # if coefficents are all 0 for too long time, shrink lambda
  tol_all0times <- ifelse(methods::hasArg(tol_all0times), dots$tol_all0times, 3)


  shrinktau2 <- ifelse(methods::hasArg(shrinktau2), dots$shrinktau2, 1.1)

  # default: no print or plot
  ifPrint <- ifelse(methods::hasArg(ifPrint), dots$ifPrint, 0)

  ########################
  # initialize
  N <- length(y) - L
  y <- y[-(1:L)]
  x <- x[-(1:L),]
  P <- D*L
  L0 <- P * spconfig$nBspline
  dIndex <- rep(1:D, L)




  lambdas <- c(1/delta_multiplier, 1, delta_multiplier) * lambda_init
  tau2 <- tau2_init

  #add optimal return result
  lambda_opt <- rep(lambda_init ,N)
  tau2_opt <- rep(tau2_init ,N)

  X0 <- matrix(NA, nrow = N, ncol = L0)
  beta_history <- matrix(NA, nrow = N, ncol = 3 * L0)
  predict1 <- matrix(NA, nrow = N, ncol = 3)
  preErr <- matrix(NA, nrow = N, ncol = 3)
  lambda_history <- matrix(NA, nrow = N, ncol = 3)

  # maintain index begin as the current beginning row index after dead loop happens
  begin <- n <- 1

  while (n <= N) {
    if (!n%%100) print(sprintf('iteration %d', n))

    if (n > (L+1) && checkDeadLoop > 5)  {
      # if dead loop happens, restart the program immediately
      warning('Dead loop happens! (runtime error) restart program ... ')
      # reset the pointer of begining row index
      begin <- n
    }

    if (n == begin) {
      ###################################
      # Initial Stage

      # initialize A and B and w
      X0[begin, ] <- getRegressor(x[begin,], dIndex, spconfig)
      reg <- X0[begin, ]
      beta <- rep(0.1 * stats::rnorm(L0), 3)
      A = forget_factor[1] * reg %*% t(reg)
      B = forget_factor[1] * reg * y[begin]

      tau2 = 0.5 / eigen(A)$values[1] # prevent EM blow up

      # record beta
      beta_history[begin, ] <- beta

      # initialize control parameter
      prepareTime <- 1 * testSize # count time (decreasing), get ready for next comparison when <= 0
      n_move_tau2 <- 1
      checkDeadLoop <- 0
      n_move_lambda <- 1
      n_last_rejuvenate <- begin
      n_rejuvenate <- n + 1 # the index for rejuvenation
      A_rejuvenate <- A # TODO: seems to me problematic because it does not match with n_rejuvenate
      B_rejuvenate <- B

    } else {
      ###################################
      # Usual Stage
      prepareTime = prepareTime - 1

      #compute new reg and update the mean of X, Y
      X0[n, ] = getRegressor(x[n,], dIndex, spconfig)
      mX = colMeans(X0[1:n,], na.rm = TRUE)
      mY = mean(y[1:n])

      reg = X0[n,] - mX

      # compute the prediction error
      for (c in 1:3) {
        predict1[n,c] = mY + sum(reg * beta[((c-1)*L0+1):(c*L0)])
        preErr[n,c] = (y[n] - predict1[n,c])^2
      }

      # update sufficient statistics
      A = (1-forget_factor[n]) * A + forget_factor[n] * reg %*% t(reg)
      B = (1-forget_factor[n]) * B + forget_factor[n] * (y[n]-mY) * reg

      if(sum(is.na(A)) != 0 || sum(is.na(B))!=0){
        print('A and B contains NA, stop and n ' )

      }

      success = rep(FALSE, 3)
      beta_old <- beta
      for (c in 1:3) {
        out <- glasso_EM(beta_old[((c-1)*L0+1):(c*L0)], A, B, D*L, lambdas[c], tau2, ...)
        success[c] <- out$success
        beta_history[n, ((c-1)*L0+1):(c*L0)] <- beta[((c-1)*L0+1):(c*L0)] <- out$beta
      }

      # force lambdas to reduce if the smallest lambda still produce all zero
      # estimates, begin after large n in order to avoid initialization issue
      if (n - begin > testSize && sum(sum(abs(beta_history[(n-tol_all0times):(n-1), 1:L0]))) == 0) {
        # if the smallest penalty produces all zero
        current_safeShrink_lambda = 1 + shrinkStepSize[n_move_lambda]/shrinkStepSize[1] * (safeShrink_lambda-1)
        lambdas = lambdas / current_safeShrink_lambda
        n_move_lambda <- n_move_lambda + 1
        if (ifPrint) {
          print(sprintf('shrink lambda(2) to %.3f because of all zero w in channel 1 at n=%d', lambdas[2], n))
        }
      }

      # ensure that tau2 is not too large such that errors blow up
      if ( sum(success) < 3 ) {
        # if blow up detected in any channel
        old_tau2 <- tau2
        tau2 = 0.5 / ( eigen(A)$values[1] * shrinktau2^(checkDeadLoop) )

        if (ifPrint) {
          print(sprintf('shrink tau2 from %.3f to %.3f with checkDeadLoop %d and n=%d rejuvenates to %d',
                        old_tau2, tau2, checkDeadLoop, n, n_rejuvenate))
        }

        n_move_tau2 = n

        # rejuvenate
        n = n_rejuvenate
        A = A_rejuvenate
        B = B_rejuvenate
        prepareTime = testSize

        # check dead loop
        checkDeadLoop <- ifelse(n_rejuvenate == n_last_rejuvenate, checkDeadLoop + 1, 0)
        n_last_rejuvenate = n_rejuvenate

        next
      }

      if (n - n_move_tau2 > tol_idle_tau2) {
        # if (ifPrint) {
        #   print(sprintf('adjust tau2 from %.4f to %.4f at n=%d', tau2, 0.5 / eigen(A)$values[1], n))
        # }
        tau2 = 0.5 / eigen(A)$values[1]
        n_move_tau2 = n
      }

      if (prepareTime <= 0) {
        # update rejuvenate point
        n_rejuvenate = n
        A_rejuvenate = A
        B_rejuvenate = B

        # switch middle channel after we have roughly accurate prediction error to determine best channel
        perfor = colMeans(preErr[(n-testSize+1):n, ])
        # get the optimal channel
        # favor larger lambda within spaTol_c tolerance
        idx_c = which.min(perfor * c(sparse_lambda_tolerance^2, sparse_lambda_tolerance, 1))

        if (!any(as.logical(beta_history[n,]))) {
          idx_c = 2
          print(sprintf('Do not move because of all zero w at n=%d',n))
        }

        if (idx_c != 2) {
          # move channel as the middle one is not the best
          if (ifPrint) {
            print(sprintf('move from c = 2(%.3f) to %d(%.3f) at %dth point with error from %.2e to %.2e',
                          lambdas[2], idx_c, lambdas[idx_c], n, perfor[2], perfor[idx_c]))
          }
        }

        # shrink the move step of lambda according to the prescribed adaptivity
        currentdelta_multiplier = 1 + shrinkStepSize[n_move_lambda]/shrinkStepSize[1] * (delta_multiplier-1)
        lambdas = lambdas[idx_c] *c(1/currentdelta_multiplier, 1, currentdelta_multiplier)
        n_move_lambda <- n_move_lambda + 1
        prepareTime = testSize


        # TODO: check if useful, may just delete them if not needed
        # if (idx_c < 2) {
        #   # it is very weird because w_path(1) would contain wrong info if reused
        #   # wish R and Matlab has pointer
        #   w_path(3) = w_path(2)
        #   w_path(2) = w_path(1)
        #   preErr(3) = preErr(2)
        #   preErr(2) = preErr(1)
        # } else if (idx_c > 2) {
        #   w_path(1) = w_path(2)
        #   w_path(2) = w_path(3)
        #   preErr(1) = preErr(2)
        #   preErr(2) = preErr(3)
        # }
      }


    }
    lambda_opt[n] <- lambdas[2]
    tau2_opt[n] <- tau2
    lambda_history[n, ] <- lambdas
    beta_opt <- beta_history[n,(L0+1):(2*L0)]
    n <- n + 1
  }

  print('==========getSequentialNonlinearModel finished=============')
  # TODO: currently I do not output optimal information as in Matlab code. Done
  result = list(beta=beta, mY=mY, beta_history = beta_history,
                preErr = preErr, predict1 = predict1,
                lambda_history = lambda_history,
                lambda_opt = lambda_opt,
                tau2_opt = tau2_opt,
                beta_opt = beta_opt,
                spconfig = spconfig,
                A = A,B = B)
  saveRDS(result,"result.RDS")
  return(result)
}
