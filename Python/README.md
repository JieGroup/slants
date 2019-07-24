# The 'slants' Python API

This package provides a Python implementation as well as an user-friendly API for the **Sequential Adaptive Nonlinear Modeling of Time Series** (SLANTS) algorithm, as proposed by **Q. Han, J. Ding, E. M. Airoldi and V. Tarokh** in the _IEEE Transactions on Signal Processing Journal, Vol. 65, NO. 19, October 2017, pages [4994, 5005]_. SLANTS provides a new method for online modeling and prediction of nonlinear and nonparametric autoregressive time series. 

- It uses **splines** to approximate a wide range of **nonlinear functions and adaptive filtering** to accommodate time varying data generating processes. 
- It is built on a **new online group LASSO algorithm** proposed in the reference paper. 
- It can be **applied to high dimensional time series** where the dimension is larger than the sample size.

Further details and theory about the algorithm can be found here [[PDF](http://jding.org/jie-uploads/2018/11/slant.pdf)]

## A Quick Setup Guide

### Getting Started 

#### 1. Install the 'slants' package using pip

```bash
# Installing test package
python -m pip install slants
```

#### 2. Import the Model and Experiment API classes

Two fundamental API classes are available. The **Model Class API** provides online model fitting capabilities of the SLANTS algorithm for your data as well as diagnostic plot graphing functions to view learned optimal model parameters. The **Experiment Class API** is a utility data-loader class that will quickly set you up with the environment to test the SLANTS algorithm API by loading sample synthetic experiments. You may import both classes as follows:

```python
from slants import Model, Experiment
```

### Using the package

#### Experiment Class API

##### 1. Load a synthetic experiment dataset using the Experiment.load() method 
This method returns a dictionary containing the design data matrix + tunable parameters for the model.

```python
exp = Experiment.load(1) # load synthetic experiment no. (1)
```
##### 2. View number of available synthetic datasets

```python
Experiment.load('?') 
```

All the tunable parameters for model fitting, inclusive of the data design matrix, can be accessed from the 'exp' dictionary.

#### Model Class API

##### 1. Instantiate Model objects and viewing available methods:

```python
model = Model()
model.help(ls=True) #ls ~ inspiration from linux bash-style command
```
##### 2. Fit SLANTS algorithm on data design matrix:

```python
# Consolidating spline spline configuration + optional tunable parameters 
spline_info = {
    'order': spline_info_order, 
    'num_Bsplines': spline_info_num_Bsplines
}

optional_params = {
    'test_size': test_size, 
    'safe_shrink_lambda': safe_shrink_lambda, 
    'random_beta': random_beta
}
## Note: Many more optional tunable parameters are available; see documentation

# Fitting the SLANTS model on the data 
model.fit(data, predict, maximum_lag, spline_info, forget_factor,  
          lambda_init, tau2_init, sparse_lambda_tolerance,		 				  
          delta_multiplier, shrink_step_size, **optional_params) # unpacking optional args dictionary 
```

Either check the **Complete_User_Guide.IPYNB** notebook or type the following to get additional documentation regarding the Model.fit() method:

```python
model.help('fit')
```

##### 3. View diagnostic plots

```python
model.plot('coef') # show trend of beta coefficients with time
model.plot('lambda') # show trend of LASSO penalty term with time
model.plot('tau') # show trend of EM decomposition parameter with time
model.plot('trace') # show trace of nonlinear functions
``` 

##### 4. Return learned model parameters

```python
learned_model_params = model.model_params() # dictionary consisting of all learned parameters
```

Either check **Complete_User_Guide.IPYNB** notebook or type the following to get additional documentation regarding the Model.model_params() method:

```python
model.help('model_params')
```

##### 5. Reset the model in rare cases of unexpected behaviour

```python
model.reset_model()
``` 

##### 6. View man-style documentation of the Model class/methods

```python
model.help(ls=True) # show list of methods in Model class 
model.help('fit') # show info about function name, for eg. fit()
model.help('all') # show info about entire class
```

## Software Guide

Interactive hands-on guides pertaining to the package may be accessed [[online](https://colab.research.google.com/drive/131JrKNapbWYMLUk67N2UUSBcFiiY_1kd)] on Google Colab or is available for download in the ```docs/``` folder of this project. Developers can avail all software-design charts of this project in the ```docs/``` folder. 

## Relevant Papers

Q .Han, J. Ding, E. Airoldi, V. Tarokh, "SLANTS: Sequential adaptive nonlinear modeling of time series," IEEE Transactions on Signal Processing 65 (19), 4994-5005. [[PDF](http://jding.org/jie-uploads/2018/11/slant.pdf)]

X. Xian, J. Ding, "Physics-assisted online learning," preprint.

## Acknowledgement

This research is funded by the Defense Advanced Research Projects Agency (DARPA) under grant number HR00111890040.

## License

The software is subjected to the GNU GPLv3 licensing terms and agreements. 