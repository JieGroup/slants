# The 'slants' package

This project repository provides software for the algorithm of Sequential Adaptive Nonlinear Modeling of Time Series developed by Jiaqi Liu and Rachit Jas. The work is an extension of a prior work of **Q. Han, J. Ding, E. M. Airoldi and V. Tarokh** in the IEEE Transactions on Signal Processing Journal, Vol. 65, NO. 19, October 2017. The open source project extends in several technical details to enhance prediction accuracy and computational efficiency. SLANTS provides a new method for online modeling and prediction of nonlinear and nonparametric autoregressive time series. 

## Available Implementations

At this time, the 'slants' API is available right out of the box to the general public for personal use in the following programming languages:
- [R](R/) - [User Guide](https://github.com/JieGroup/slants/blob/master/R/vignettes/user-guide.pdf)
- [Python](https://pypi.org/project/slants/) - [User Guide](https://colab.research.google.com/drive/131JrKNapbWYMLUk67N2UUSBcFiiY_1kd)

## Main Features of the Algorithm

- It uses **splines** to approximate a wide range of **nonlinear functions and adaptive filtering** to accommodate time varying data generating processes. 
- It is built on a **new online group LASSO algorithm** proposed in the reference paper. 
- It can be **applied to high dimensional time series** where the dimension is larger than the sample size.

## Relevant Papers

Q .Han, J. Ding, E. Airoldi, V. Tarokh, "SLANTS: Sequential adaptive nonlinear modeling of time series," IEEE Transactions on Signal Processing 65 (19), 4994-5005. [[PDF](http://jding.org/jie-uploads/2018/11/slant.pdf)]

J. Liu, J. Zhou, J. Ding, "Privacy aware supervised learning," preprint.

X. Xian, J. Ding, "Physics-assisted online learning," preprint.



## License

The software is subjected to the GNU GPLv3 licensing terms and agreements. 
