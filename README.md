# The 'slants' R package
This package implements and extends the algorithm  'Sequential Adaptive Nonlinear Modeling of Time Series (SLANTS)'. 
SLANTS is a new method for online modeling and prediction of nonlinear and nonparametric autoregressive time series. It uses splines to approximate a wide range of nonlinear functions and adaptive filtering to accommodate time varying data generating processes. It is built on a new online group LASSO algorithm proposed in the reference paper. It can be applied to high dimensional time series where the dimension is larger than the sample size. 

## Getting Started

First install the devtools package

install.packages("devtools")

library("devtools")

Then install this package

install_github('JieGroup/slants')

## Using This Package

To see the available function to use, type 

ls("package:slants")

A quick guide of package can be found [here](https://github.com/JieGroup/slants/blob/master/vignettes/user-guide.pdf) 


## Relevant Papers

Q .Han, J. Ding, E. Airoldi, V. Tarokh, "SLANTS: Sequential adaptive nonlinear modeling of time series," IEEE Transactions on Signal Processing 65 (19), 4994-5005. [pdf](http://jding.org/jie-uploads/2018/11/slant.pdf)

X. Xian, J. Ding, "Physics-assisted online learning," preprint.

## Acknowledgment

This research is funded by the Defense Advanced Research Projects Agency (DARPA) under grant number HR00111890040.


