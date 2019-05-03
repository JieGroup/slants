pkgname <- "slants"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "slants-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('slants')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("augknt")
### * augknt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: augknt
### Title: augknt
### Aliases: augknt

### ** Examples

getRegressor()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("augknt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getPreprocess")
### * getPreprocess

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getPreprocess
### Title: get Preprocess
### Aliases: getPreprocess

### ** Examples

getPreprocess()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getPreprocess", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getRegressor")
### * getRegressor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getRegressor
### Title: get Regressor
### Aliases: getRegressor

### ** Examples

x =  matrix(rnorm(2*(1000+8)),1000+8, 2)
knotBox = apply(x, 2, quantile, c(0.02,0.98), na.rm = TRUE)
knots = apply(knotBox, 2, augknt, nBspline, order)
spconfig = list(order = 3, nBspline = 10, knots = knots,knotBox = knotBox)
regressor <- getRegressor(x[,1],1,spconfig)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getRegressor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getSequentialNonlinearModel")
### * getSequentialNonlinearModel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getSequentialNonlinearModel
### Title: getSequentialNonlinearModel
### Aliases: getSequentialNonlinearModel
### Keywords: linear model non sequential series time

### ** Examples

utils.lag <- function(ts, lag = 1, pad = NA) {
# return the lagged version of a time series vector
return(c(rep(pad, lag), ts[1:(length(ts)-lag)]))
}
# pack configurations so that it is easier to bundle things together
Ex1 <- (function(N=2000, D=2, L=8, ifPlot = TRUE){
 err <- matrix(rnorm(D*(N+L)), N+L, D)
 X <- err[ ,1]
 X <- cbind(X, 0.5 * utils.lag(X,1)^2 - 0.8 * utils.lag(X,7) + 0.2 * err[,2])
 if (ifPlot) {
   plot(X[,1], X[,2], pch = 16, cex = 0.5, col = "red")
 }
 list(N=N, D=D, L=L, X=X[-(1:L),], y=X[-(1:L),2])
})()

Ex1_algo <- (function(experiment_config){
ec <- experiment_config # short name
lambda <- 1/c(1:ec$N) # same as batch
shrinkStepSize <- 1/c(1:ec$N)
list(
  order = 3,
  nBspline = 10,
  lambda = lambda,
  shrinkStepSize = shrinkStepSize,
  # if the performance does not exceed this ratio, the sparser (larger gamma) is preferred
  spaTol_gamma = 1.01,
  moveSize = 10^(0.4), # multipler to move among channels
 safeShrink_gamma = 10^(0.1), # to adjust gamma before they get too crazy
  gamma_init = 0.01,
  alpha2_init = 0.05
)
})(Ex1)

Ex1_algo <- c(Ex1_algo, do.call(getPreprocess, c(Ex1, Ex1_algo)))

Ex1_result <- do.call(getSequentialNonlinearModel, c(list(ifPrint=TRUE, testSize=50),Ex1, Ex1_algo))
========diagonistic========
plot(Ex1_result$gamma_opt,type = "l")
plot(Ex1_result$alpha_opt,type = "l",ylab = "Tao2")
plot(Ex1_result$preErr[,2],type = "l")
#get the historical opt beta, always the middle channel
beta_hist_opt <- Ex1_result$beta_history[,161:320]

order <- 3
knots <- Ex1_algo$knots[,1]
for (t in length(beta_hist)){
tmp <- NULL
for (i in 1:(length(knots)-order)) {
 tmp <- cbind(tmp, sapply(seq(-1,1,by=0.01), spval, x=knots[i:(i+order)]))
}
par(mfrow = c(1,2))
plot(seq(-1,1,by=0.01), tmp %*% beta_hist_opt[t,1:10])
plot(seq(-1,1,by=0.01), tmp %*% beta_hist_opt[t,121:130])
par(mfrow=c(1,1))
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getSequentialNonlinearModel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("glasso_EM")
### * glasso_EM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glasso_EM
### Title: glasso_EM
### Aliases: glasso_EM

### ** Examples

getRegressor()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glasso_EM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotcoeff")
### * plotcoeff

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotcoeff
### Title: plotcoeff
### Aliases: plotcoeff

### ** Examples

positive.polynomial(0.2,x = 0,k = 3,end_knot = 2,d=0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotcoeff", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("positive.polynomial")
### * positive.polynomial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: positive.polynomial
### Title: positive.polynomial
### Aliases: positive.polynomial

### ** Examples

positive.polynomial(0.2,x = 0,k = 3,end_knot = 2,d=0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("positive.polynomial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spval")
### * spval

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spval
### Title: divided difference
### Aliases: spval

### ** Examples

getRegressor()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spval", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
