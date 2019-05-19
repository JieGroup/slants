#' getSequentialNonlinearModel
#'
#' This function is to fit seqential non linear model
#'
#' @param y vector
#' @param x the reuslt from getPreprocess \code{x}
#' @param D columns of original X
#' @param L maximize lag number
#' @param lambda the LASSO penalty
#' @param gamma_init the initial value for gamma
#' @param alpha2_init the initial value of EM decomposition parameter, refer to the paper
#' @param spconfig spline configuration
#' @param moveSize multipler to move among channels
#' @param spaTol_gamma favor smaller gamma within spaTol_c tolerance
#' @param shrinkStepSize adjust gamma before they get too crazy
#' @param ... default parameters,including
#'  \itemize{
#'   \item testSize: the number of obs to average prediction error, default is 50
#'   \item tol_idle_alpha2:increase alpha2 if it idles for a while
#'   \item safeShrink_gamma : adjust gamma before they get too crazy
#'   \item tol_all0times: if coefficents are all 0 for too long time, shrink gamma
#'   \item shrinkAlpha2: how to shrink alpha2
#' }
#' @return a list including
#'  \itemize{
#'   \item preErr: historical predicted error, firt column indicates the first channel, etc.
#'   \item predict1 : historical predicted values, firt column indicates the first channel, etc.
#'   \item gamma_history : historical gamma channels, first channel is the smallest one and the third channel is the largest one.
#'   \item gamma_opt: optimized gamma, the second channel in gamma_history
#'   \item alpha_opt: optimized alpha
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
#' lambda <- 1/c(1:ec$N) # same as batch
#' shrinkStepSize <- 1/c(1:ec$N)
#' list(
#'   order = 3,
#'   nBspline = 10,
#'   lambda = lambda,
#'   shrinkStepSize = shrinkStepSize,
#'   # if the performance does not exceed this ratio, the sparser (larger gamma) is preferred
#'   spaTol_gamma = 1.01,
#'   moveSize = 10^(0.4), # multipler to move among channels
#'  safeShrink_gamma = 10^(0.1), # to adjust gamma before they get too crazy
#'   gamma_init = 0.01,
#'   alpha2_init = 0.05
#' )
#' })(Ex1)
#'
#'Ex1_algo <- c(Ex1_algo, do.call(getPreprocess, c(Ex1,Ex1_algo)))
#'par <- c(list(ifPrint=1,testSize=50),Ex1,Ex1_algo)
#'Ex1_result <- do.call(getSequentialNonlinearModel, par)
#'#========diagonistic========
#'graphics::plot(Ex1_result$gamma_opt,type = "l")
#'graphics::plot(Ex1_result$alpha_opt,type = "l",ylab = "Tao2")
#'graphics::plot(Ex1_result$preErr[,2],type = "l")
#'#get the historical opt beta, always the middle channel
#'plotcoeff(Ex1_result$beta_opt,Ex1_algo$knots,Ex1_algo$nBspline)



getSequentialNonlinearModel <- function(y, x, D, L, lambda, gamma_init, alpha2_init, spconfig, spaTol_gamma, shrinkStepSize, moveSize, ...) {



  ########################
  # set default parameters


  dots <- list(...)

  spaTol_gamma <- ifelse(methods::hasArg(spaTol_gamma), spaTol_gamma, 1.01)

  gamma_init <- ifelse(methods::hasArg(gamma_init), gamma_init, 0.03)

  alpha2_init <- ifelse(methods::hasArg(alpha2_init), alpha2_init, 0.05)

  # testSize is the number of obs to average prediction error, default is 50
  testSize <- ifelse(methods::hasArg(testSize), dots$testSize, 50)

  # increase alpha2 if it idles for a while
  tol_idle_alpha2 <- ifelse(methods::hasArg(tol_idle_alpha2), dots$tol_idle_alpha2, 50)

  # to adjust gamma before they get too crazy
  safeShrink_gamma <- ifelse(methods::hasArg(safeShrink_gamma), dots$safeShrink_gamma, 10^(0.4))

  # if coefficents are all 0 for too long time, shrink gamma
  tol_all0times <- ifelse(methods::hasArg(tol_all0times), dots$tol_all0times, 3)


  shrinkAlpha2 <- ifelse(methods::hasArg(shrinkAlpha2), dots$shrinkAlpha2, 1.1)

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




  gammas <- c(1/moveSize, 1, moveSize) * gamma_init
  alpha2 <- alpha2_init

  #add optimal return result
  gamma_opt <- rep(gamma_init ,N)
  alpha2_opt <- rep(alpha2_init ,N)

  X0 <- matrix(NA, nrow = N, ncol = L0)
  beta_history <- matrix(NA, nrow = N, ncol = 3 * L0)
  predict1 <- matrix(NA, nrow = N, ncol = 3)
  preErr <- matrix(NA, nrow = N, ncol = 3)
  gamma_history <- matrix(NA, nrow = N, ncol = 3)

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
      A = lambda[1] * reg %*% t(reg)
      B = lambda[1] * reg * y[begin]

      alpha2 = 0.5 / eigen(A)$values[1] # prevent EM blow up

      # record beta
      beta_history[begin, ] <- beta

      # initialize control parameter
      prepareTime <- 1 * testSize # count time (decreasing), get ready for next comparison when <= 0
      n_move_alpha2 <- 1
      checkDeadLoop <- 0
      n_move_gamma <- 1
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
      A = (1-lambda[n]) * A + lambda[n] * reg %*% t(reg)
      B = (1-lambda[n]) * B + lambda[n] * (y[n]-mY) * reg

      if(sum(is.na(A)) != 0 || sum(is.na(B))!=0){
        print('A and B contains NA, stop and n ' )

      }

      success = rep(FALSE, 3)
      beta_old <- beta
      for (c in 1:3) {
        out <- glasso_EM(beta_old[((c-1)*L0+1):(c*L0)], A, B, D*L, gammas[c], alpha2, ...)
        success[c] <- out$success
        beta_history[n, ((c-1)*L0+1):(c*L0)] <- beta[((c-1)*L0+1):(c*L0)] <- out$beta
      }

      # force gammas to reduce if the smallest gamma still produce all zero
      # estimates, begin after large n in order to avoid initialization issue
      if (n - begin > testSize && sum(sum(abs(beta_history[(n-tol_all0times):(n-1), 1:L0]))) == 0) {
        # if the smallest penalty produces all zero
        current_safeShrink_gamma = 1 + shrinkStepSize[n_move_gamma]/shrinkStepSize[1] * (safeShrink_gamma-1)
        gammas = gammas / current_safeShrink_gamma
        n_move_gamma <- n_move_gamma + 1
        if (ifPrint) {
          print(sprintf('shrink gamma(2) to %.3f because of all zero w in channel 1 at n=%d', gammas[2], n))
        }
      }

      # ensure that alpha is not too large such that errors blow up
      if ( sum(success) < 3 ) {
        # if blow up detected in any channel
        old_alpha2 <- alpha2
        alpha2 = 0.5 / ( eigen(A)$values[1] * shrinkAlpha2^(checkDeadLoop) )

        if (ifPrint) {
          print(sprintf('shrink alpha2 from %.3f to %.3f with checkDeadLoop %d and n=%d rejuvenates to %d',
                        old_alpha2, alpha2, checkDeadLoop, n, n_rejuvenate))
        }

        n_move_alpha2 = n

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

      if (n - n_move_alpha2 > tol_idle_alpha2) {
        # if (ifPrint) {
        #   print(sprintf('adjust alpha2 from %.4f to %.4f at n=%d', alpha2, 0.5 / eigen(A)$values[1], n))
        # }
        alpha2 = 0.5 / eigen(A)$values[1]
        n_move_alpha2 = n
      }

      if (prepareTime <= 0) {
        # update rejuvenate point
        n_rejuvenate = n
        A_rejuvenate = A
        B_rejuvenate = B

        # switch middle channel after we have roughly accurate prediction error to determine best channel
        perfor = colMeans(preErr[(n-testSize+1):n, ])
        # get the optimal channel
        # favor larger gamma within spaTol_c tolerance
        idx_c = which.min(perfor * c(spaTol_gamma^2, spaTol_gamma, 1))

        if (!any(as.logical(beta_history[n,]))) {
          idx_c = 2
          print(sprintf('Do not move because of all zero w at n=%d',n))
        }

        if (idx_c != 2) {
          # move channel as the middle one is not the best
          if (ifPrint) {
            print(sprintf('move from c = 2(%.3f) to %d(%.3f) at %dth point with error from %.2e to %.2e',
                          gammas[2], idx_c, gammas[idx_c], n, perfor[2], perfor[idx_c]))
          }
        }

        # shrink the move step of gamma according to the prescribed adaptivity
        currentMoveSize = 1 + shrinkStepSize[n_move_gamma]/shrinkStepSize[1] * (moveSize-1)
        gammas = gammas[idx_c] *c(1/currentMoveSize, 1, currentMoveSize)
        n_move_gamma <- n_move_gamma + 1
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
    gamma_opt[n] <- gammas[2]
    alpha2_opt[n] <- alpha2
    gamma_history[n, ] <- gammas
    beta_opt <- beta_history[n,(L0+1):(2*L0)]
    n <- n + 1
  }

  print('==========getSequentialNonlinearModel finished=============')
  # TODO: currently I do not output optimal information as in Matlab code. Done
  list(mY=mY, beta_history = beta_history,
       preErr = preErr, predict1 = predict1,
       gamma_history = gamma_history,
       gamma_opt = gamma_opt,
       alpha_opt = alpha2_opt,
       beta_opt = beta_opt)
}
