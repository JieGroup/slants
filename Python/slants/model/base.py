"""SLANTS ALGORITHM IMPLEMENTATION
Utility class for data preprocessing, visualization and model fitting methods    

    Copyright (C) 2019 - SLANTS authors

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
"""

# TODO: Might need to handle return False from all functions + add variadic tuple args + abstract func args handling code


# Importing Dependencies
import pickle
import numpy as np
import warnings
from copy import deepcopy 
from pandas import DataFrame
from os.path import exists
from slants.model.spline import BSpline 
import matplotlib.pyplot as plt


class SLANTS:
	"""Utility Class for the Raw Implementation of the SLANTS algorithm

	Description:
		Provisions the underlying logic behind the building blocks of the algorithm, such as the data preprocessing and Glasso EM steps  
	""" 	 
	@staticmethod 
	def get_regressor(**args):
		"""Preprocesses original data to get splines configuration and input data for B-Splines

		Expected Args:
			x: Input data vector x
			d: Vector of same length as x, holding the index of which knots to look at
			spconfig: Splines configuration

		Returns:
			sp: Vector of regressors that will be passed to EM algorithm

		"""  	
		# Handling missing arguments
		required_args = ['x', 'd', 'spconfig']
		if set(required_args).intersection(set(args)) != set(required_args):  		
			print('[ERROR] Expected 3 required arguments in SLANTS.get_regressor(): \
				\n\t(1) x = Input data vector x \
				\n\t(2) d = Vector of same length as x, holding the index of which knots to look at\
				\n\t(3) spconfig = Splines configuration')
			return False

		# Storing required arguments
		x, d, spconfig = args['x'], args['d'], args['spconfig']
		knots, order, nBspline = spconfig['knots'], spconfig['order'], spconfig['nBspline']

		# Computing divided differenced scalars, i.e., regressors obtained by computing spline basis-function values of the input
		# Regressors computed from each input data 'x[i]' of an input vector at a time   
		sp = np.array([BSpline.spval(t=x[i], x=knots[j:(j + order) + 1, d[i]]) \
			for i in range(len(x)) for j in range(nBspline)])

		return sp    		 



	@staticmethod
	def preprocess_data(**args):
		"""Transform original data of dimension (N*D) to (N*DL) + create spline configuration for the new data

		Expected Args:
			X: Original input data of shape (N*D)
			L: Maximum lag
			order: Order of B-Spline
			nBspline: Number of B-Spline 
			augkntFlag: Binary indictor to highlight whether to find knots 
			scaleFlag: Binary indictor to highlight whether to standardize the columns of input data, X

		Returns:
			A Dictionary containing the following items:-
			x: Transformed data of shape (N*DL)
			scale: Standard deviation used for standardization of input data, X
			feasibleBox: Region of input data, X where the spline fit is very accurate
			knotBox: Knot range in the 'rescale scale', which is slightly larger than feasibleBox in quantile
			knots: Knots of the spline
			spconfig: Spline configuration for B-Splines

		"""
		# Handling missing arguments
		required_args = ['X', 'L', 'order', 'nBspline']
		if set(required_args).intersection(set(args)) != set(required_args):
			print('[ERROR] Expected 4 required arguments in SLANTS.preprocess_data(): \
				\n\t(1) X = Originalinput data of shape (N*D) \
				\n\t(2) L = Maximum Lag\
				\n\t(3) order = Order of B-Spline\
				\n\t(4) nBspline = Number of B-Spline')
			return False
		scaleFlag = True if 'scaleFlag' not in args else args['scaleFlag']
		augkntFlag = True if 'augkntFlag' not in args else args['augkntFlag']

		# Storing required arguments
		X, L, order, nBspline = args['X'], args['L'], args['order'], args['nBspline'] 

		# Storing original dimension of input data
		N, D = X.shape

		# Declaring transformed input data matrix
		x = np.zeros((D*L, N))
		x.fill(np.nan)

		# Computing standardization scaling factor for each column, if applicable  	 
		scale = np.std(X, axis=0, ddof=1, keepdims=True) if scaleFlag else np.repeat(1, D).reshape(1, -1)

		# Standardizing and arranging the data matrix in row order: X(1,n-1)...X(D,n-1), ..., X(1,n-L)...X(D,n-L) 
		for lag in range(L): x[((lag)*D):((lag+1)*D), (L):(N+1)] = X[(L-lag-1):(N-lag-1), ].T/scale.T
		x = x.T

		# Determining the feasible box in original scale - the region of X where the spline fit is very accurate
		feasibleBox = np.vstack((np.nanpercentile(X, 10, axis=0),np.nanpercentile(X, 90, axis=0)))

		# Determining rhe knot range in rescaled scale, which is slightly larger than feasible box in quantile
		knotBox = np.vstack((np.nanpercentile(x[:,:D], 2, axis=0),np.nanpercentile(x[:,:D], 98, axis=0)))

		# Determining the knots and spline basis
		if augkntFlag:
			#knots = np.apply_along_axis(BSpline.augknt, 0, knotBox, nBspline, order).reshape(-1, knotBox.shape[1])
			knots = np.apply_along_axis(lambda boundary: BSpline.augknt(boundary=boundary, nBspline=nBspline, order=order), axis=0, arr=knotBox).reshape(-1, knotBox.shape[1])
		else:
			knots = knotBox[[0, 1, 1],:]

		# Structuring the spline configuration
		spconfig = {
			'order': order,
			'nBspline': nBspline,
			'knots': knots,
			'knotBox': knotBox
		}

		# Returns the transformed data + sample standard deviation + spline configuration
		return dict(x=x, scale=scale.reshape(-1), feasibleBox=feasibleBox, knotBox=knotBox, knots=knots, spconfig=spconfig) #Note: scale is a 2D np.array



	@staticmethod
	def fit(**args): #getSequentialNonlinearModel.R
		"""Fits Sequential Non Linear Model to Time Series data

		Expected Args:
			y: Original input data of shape (N*D)
			x: Preprocessed input data, x
			D: Columns of original data, X
			L: Maximum Lag
			lambda_value: forgetting factor
			gamma_init: Initial value for lasso penalty  
			alpha2_init: Initial value for EM decomposition parameter, tao^2
			spconfig: Spline configuration 
			move_size: Multiplier to move among channels 
			spaTol_gamma: Favors smaller gamma within spaTol_c tolerance 
			shrink_step_size: Adjusts gamma before they get too crazy
			y_index: Regressor in consideration 
		
		Default Args:
			test_size: Number of observations to average prediction error, default 50 
			tol_idle_alpha2: Increase alpha2 if it idles for a while 
			safe_shrink_gamma: Adjust gamma before thy get too crazy 
			tol_all0times: If coefficients = 0 for a long time, then shrink gamma 
			shrink_alpha2: If alpha2 is too big, then shrink it 
			debug_print: Don't print or plot by default, 0
			use_random_beta: Provide betas for initialization
			online: Use online version of sequential model using serialized model params, default True.  

		Returns:
			A Dictionary containing the following items:-
			pre_err: historical predicted error; 1st column indicates 1st channel, etc. 
			predict1: historical predicted values; 1st column indicates 1st channel, etc. 
			gamma_history: historical gamma channels; 1st channel is the smallest one, 3rd channel is the largest one   
			gamma_opt: optimized gamma; 2nd channel in gamma_history 
			alpha_opt: optimized alpha 
			beta_opt: best coefficient, the last row and second channel related in beta_history  
			beta_history: all stored coefficient
			spconfig: spline configuration
			A, B: sufficient statistics

		"""
		# Handling missing arguments
		required_args = ['y', 'x', 'D', 'L', 'lambda_value', 'gamma_init', 'alpha2_init', 'spconfig', 'spaTol_gamma', 'shrink_step_size', 'move_size']
		if set(required_args).intersection(set(args)) != set(required_args):
			print('[ERROR] Expected 11 required arguments in SLANTS.fit(): \
				\n\t(1) y = Original input data of shape (N*D)\
				\n\t(2) x = Preprocessed input data, x \
				\n\t(3) D = Columns of original data, X \
				\n\t(4) L = Maximum Lag \
				\n\t(5) lambda_value = LASSO penalty \
				\n\t(6) gamma_init = Initial value for gamma \
				\n\t(7) alpha2_init = Initial value for EM decomposition parameter, tao^2 \
				\n\t(8) spconfig = Spline configuration \
				\n\t(9) move_size = Multiplier to move among channels \
				\n\t(10) spaTol_gamma = Favors smaller gamma within spaTol_c tolerance \
				\n\t(11) shrink_step_size = Adjusts gamma before they get too crazy ')
			return False
		y, x, D, L, lambda_value, gamma_init, alpha2_init, spconfig, spaTol_gamma, shrink_step_size, move_size = (args[element] for element in required_args)		
		test_size = 50 if 'test_size' not in args else args['test_size']
		tol_idle_alpha2 = 50 if 'tol_idle_alpha2' not in args else args['tol_idle_alpha2']
		safe_shrink_gamma = 10**(0.4) if 'safe_shrink_gamma' not in args else args['safe_shrink_gamma']
		tol_all0times = 3 if 'tol_all0times' not in args else args['tol_all0times']
		shrink_alpha2 = 1.1 if 'shrink_alpha2' not in args else args['shrink_alpha2']
		debug_print = False if 'debug_print' not in args else args['debug_print']
		use_random_beta = False if 'use_random_beta' not in args else args['use_random_beta']
		online = True if 'online' not in args else args['online'] #/////////////////////
		y_index = "_" if 'y_index' not in args else args['y_index']
		warnings.filterwarnings("ignore", category=np.ComplexWarning)

		######################## INITIALIZATION STAGE
		N = len(y) - L # new size after excluding 1st L rows
		y = np.delete(y, (np.round(np.linspace(1,L,L).astype(int)) - 1).tolist(), axis=0) # Exclude 1:L rows, first L rows
		x = np.delete(x, (np.round(np.linspace(1,L,L).astype(int)) - 1).tolist(), axis=0)
		P = D*L
		L0 = P*spconfig['nBspline']
		d_index = np.array(list(np.arange(D))*L)

		if online:
			if  exists('model-params_predict_random-variable_{}.pickle'.format(y_index+1)):
				# print('reading', flush=True)
				with open('model-params_predict_random-variable_{}.pickle'.format(y_index+1), 'rb') as serde_handle:
					history = pickle.load(serde_handle)
			else:
				# print('not reading', flush=True)
				del args['online']			
				return SLANTS.fit(**args, online=False)

		if online:
			length_1 = history['beta_history'].shape[0]
			beta = np.copy(history['beta_history'][length_1-1, :])
			hist_spconfig = deepcopy(history['spconfig'])
			lenn = len(beta)
			A, B = np.copy(history['A']), np.copy(history['B'])

			if len(x[0,:]) != lenn/(hist_spconfig['nBspline']*3): print('<WARNING> Input values need to be preprocessed first')							

		gammas = np.array([1/move_size, 1, move_size]) * gamma_init # in paper, the lambdas for each channel
		alpha2 = alpha2_init

		rev_N = N+length_1 if online else N
		gamma_opt = np.repeat(gamma_init, rev_N).astype('complex')
		alpha_opt = np.repeat(alpha2_init, rev_N).astype('complex')

		X0 = np.empty((N, L0)); X0.fill(np.nan)
		beta_history = np.empty((rev_N, 3 * L0)); beta_history.fill(np.nan); 
		if online: beta_history[0:length_1,:] = np.copy(history['beta_history']) 
		predict1 = np.empty((rev_N, 3)); predict1.fill(np.nan)
		if online: predict1[0:length_1,:] = np.copy(history['predict1']) 		
		pre_err = np.empty((rev_N, 3)); pre_err.fill(np.nan)
		if online: pre_err[0:length_1,:] = np.copy(history['pre_err']) 		
		gamma_history = np.empty((rev_N, 3)).astype('complex'); gamma_history.fill(np.nan)
		if online: gamma_history[0:length_1,:] = np.copy(history['gamma_history']) 


		# Maintain index 'begin' as the current beginning row after dead loop happens
		begin = n = 1

		while n <= N:
			if not n%100: print('Iteration: {}'.format(n))

			if n > (L+1) and check_deadloop > 5:
				# If dead loop happens, restart the program immediately
				print('Dead loop happens! (runtime error) restart program ...')
				# reset the pointer of beginning row index
				begin = n

			if n == begin:
				######################## INITIAL STAGE

				if not online:
					# Initialize A,B and w
					X0[begin-1,:] = SLANTS.get_regressor(x=x[begin-1,:], d=d_index, spconfig=spconfig)
					reg = np.array([X0[begin-1,:]])
					if type(use_random_beta) == bool:
						beta = np.repeat(0.1 * np.random.randn(L0), 3)
					else:
						beta = np.copy(use_random_beta)
					A = lambda_value[0] * np.dot(reg.T, reg) # W_T is a diagonal matrix
					B = lambda_value[0] * reg * y[begin-1]


					# Prevent EM blow up
					eigen_values, _ = np.linalg.eig(A)
					alpha2 = 0.5 / eigen_values[0] #alpha2 = alpha_init?

					# Record beta
					beta_history[begin-1, :] = np.copy(beta)

				# Initialize control parameter
				# n_* refers to index parameters 
				prepare_time = 1 * test_size # count time (decreasing), get ready for next comparison when <= 0
				n_move_alpha2, n_move_gamma = 1, 1				
				n_last_rejuvenate, n_rejuvenate = begin, n + 1  # index for rejuvenation CHECKK
				check_deadloop = 0
				A_rejuvenate, B_rejuvenate = np.copy(A), np.copy(B)

				if online:
					for c in range(0, 3):
						predict1[n+length_1-1,c] = predict1[length_1-1,c]
						pre_err[n+length_1-1,c] = pre_err[length_1-1,c]

			else:
				######################## USUAL STAGE
				prepare_time = prepare_time - 1

				# Compute new regressor and update the mean of X, Y
				X0[n-1,:] = SLANTS.get_regressor(x=x[n-1, :], d=d_index, spconfig=spconfig)
				mX = np.nanmean(X0[0:n, :], axis=0)
				mY = np.mean(y[0:n])

				reg = np.array([X0[n-1,:] - mX])

				# Compute the prediction error
				for c in range(3):
					indices = list(range(((c+1)-1)*L0, (c+1)*L0))
					if online:
						predict1[n+length_1-1, c] = mY + np.sum(reg * beta[indices].T)
						pre_err[n+length_1-1, c] = (y[n-1] - predict1[n-1, c])**2
					else:
						predict1[n-1, c] = mY + np.sum(reg * beta[indices].T)
						pre_err[n-1, c] = (y[n-1] - predict1[n-1, c])**2

				# Update sufficient statistics
				A = (1-lambda_value[n-1]) * A + lambda_value[n-1] * np.dot(reg.T, reg) # 1 - gamma*A ... in paper
				B = (1-lambda_value[n-1]) * B + lambda_value[n-1] * (y[n-1]-mY) * reg

				if np.sum(np.isnan(A)) or np.sum(np.isnan(B)): # Checking if A or B have NaN elements
						print('A and/or B contains NA, stop and n')

				# Apply Glasso EM and track success
				success = np.repeat(False, 3)
				beta_old = np.copy(beta)
				for c in range(3):
					indices = list(range(((c+1)-1)*L0, (c+1)*L0))
					out = SLANTS.apply_glasso_EM(beta_init=np.array([beta_old[indices]]), A=A, B=B, ngroup=D*L, lambda_value=gammas[c], alpha2=alpha2, debug_print=debug_print)
					success[c] = deepcopy(out['success']) # status of EM algorithm run
					beta[indices] = np.copy(out['beta']) #
					beta_history[n-1, indices] = np.copy(beta[indices])

				# Force gammas to reduce if the smallest gamma still produces all zero estimates,
				# Begin after large n in order to avoid initialization issue

				if (n - begin) > test_size and not np.sum(np.abs(beta_history[(n-1 - tol_all0times)-1 : (n-1), 0:L0])):
					current_safe_shrink_gamma = 1 + shrink_step_size[n_move_gamma-1]/shrink_step_size[0] + (safe_shrink_gamma-1)
					gammas = gammas/ current_safe_shrink_gamma
					n_move_gamma += 1
					if debug_print:
						print('Shrink gamma(channel-2) to {} because of all zero w in channel 1 at n={}'.format(gamma[1], n))

				# Ensure that alpha is not too large such that errors blow up
				if np.sum(success) < 3:
					old_apha2 = alpha2
					eigen_values, _ = np.linalg.eig(A)
					alpha2 = 0.5 / eigen_values[0] / shrink_alpha2**(check_deadloop)

					if debug_print:
						print('Shrink alpha2 from {} to {} with check_deadloop {} and n={} rejuvenates to {}'.format(old_apha2, alpha2, check_deadloop, n, n_rejuvenate))

					n_move_alpha2 = n

					# Rejuvenate
					n, A, B = n_rejuvenate, np.copy(A_rejuvenate), np.copy(B_rejuvenate)
					prepare_time = test_size

					# Check dead loop
					check_deadloop = check_deadloop + 1 if n_rejuvenate == n_last_rejuvenate else 0
					n_last_rejuvenate = n_rejuvenate

					continue

				if (n - n_move_alpha2) > tol_idle_alpha2:
					eigen_values, _ = np.linalg.eig(A)
					alpha2 = 0.5 / eigen_values[0]
					n_move_alpha2 = n

				if prepare_time <= 0:
					# Update rejuvenate point
					n_rejuvenate, A_rejuvenate, B_rejuvenate = n, np.copy(A), np.copy(B)

					# Switch middle channel after we have roughly accurate prediction error to determine best channel
					perfor = np.mean(pre_err[((n-1)-test_size+1):n, :], axis=0)

					# Get the optimal channel, favour larger gamma within spaTol_c tolerance
					idx_c = np.argmin(perfor * np.array([spaTol_gamma**2, spaTol_gamma, 1]))

					if not np.any(beta_history[n-1, :]):
						idx_c = 1
						print('Do not move because all of zeros at n={}'.format(n-1))

					if idx_c != 1 and debug_print:
						print('Move from c=1({}) to {}({}) at {}th point with error from {} to {}'.format(gammas[1], idx_c, gammas[idx_c], n , perfor[1], perfor[idx_c]))

					# Shrink the move step of gamma according to the prescribed adaptivity
					current_move_size = 1 + shrink_step_size[n_move_gamma-1]/shrink_step_size[0]
					gammas = gammas[idx_c] * np.array([1/current_move_size, 1, current_move_size])
					n_move_gamma = n_move_gamma + 1
					prepare_time = test_size					

			######################## COMPLETION STAGE (out of the loop)
			if online:
				# Storing optimal parameters
				gamma_opt[lenn+n-1] = np.copy(gammas[1])
				alpha_opt[lenn+n-1] = np.copy(alpha2)
				gamma_history[lenn+n-1,:] = np.copy(gammas) 
				beta_opt = np.copy(beta_history[lenn+n-1, (L0):(2*L0)])
			else:
				# Storing optimal parameters
				gamma_opt[n-1] = np.copy(gammas[1])
				alpha_opt[n-1] = np.copy(alpha2)
				gamma_history[n-1,:] = np.copy(gammas) 
				beta_opt = np.copy(beta_history[n-1, (L0):(2*L0)])
			n += 1

		print('==========< Sequential Nonlinear Model finished >=============')
		
		# Serializing mode parameters
		model_params = dict(
							mY=mY, beta_history=beta_history, pre_err=pre_err, predict1=predict1,
							gamma_history=gamma_history, gamma_opt=gamma_opt, alpha_opt=alpha_opt,
							beta_opt=beta_opt, spconfig=spconfig, A=A, B=B
						)
		with open('model-params_predict_random-variable_{}.pickle'.format(y_index+1), 'wb') as serde_handle:
			pickle.dump(model_params, serde_handle)

		return model_params



	@staticmethod
	def plot_coeff(**args):
		"""Plots the trace of nonlinear functions

		Expected Args:
			beta: Optimal Beta coefficients
			knots: Knots returned from preprocess_data() 
			nBspline: Number of B-splines needed to fit the model
			order: Order of splines
		"""
		# Handling missing arguments
		required_args = ['beta', 'knots', 'nBspline', 'order']
		if set(required_args).intersection(set(args)) != set(required_args):
			print('[ERROR] Expected 4 required arguments in SLANTS.plot_coeff(): \
				\n\t(1) beta = Optimal Beta coefficients\
				\n\t(2) knots =  Knots, returned from preprocess_data() \
				\n\t(3) nBspline = Number of B-splines needed to fit the model \
				\n\t(4) order = Order of splines')
			return False
		beta, knots, nBspline, order = (args[element] for element in required_args)

		# Preparing primary data attributes
		knots = np.copy(knots[:,0]) # make knot 1D
		D = len(knots)
		DL = len(beta)//nBspline
		x = np.arange(-1, 1.01, 0.01)
		tmp = None

		# Computing tmp and its stacking columns vertically side by side
		for i in range(len(knots)-order):
			col = np.array([BSpline.spval(t=np.array(t), x=knots[i:(i+order)+1]) for t in np.arange(-1, 1.01, 0.01)])
			if tmp is None: tmp = col
			else: tmp = np.column_stack((tmp, col))

		# Initializing result matrix
		result = DataFrame(np.nan, index=['lag{}'.format(lag) for lag in range(1, DL+1)], columns=[str(col_index) for col_index in range(1,tmp.shape[0]+1)])
		result.fillna(method='ffill', inplace=True)

		# Populating result matrix 
		for i in range(DL):
			beta_new = np.copy(beta[(nBspline *(i)):(nBspline*(i+1))])
			result.iloc[i,:] = np.dot(tmp, beta_new)

		# Setting up figure and subplot
		fig = plt.figure()
		ax = plt.subplot(111)

		# Plotting figures
		ax.plot(x, result.iloc[0,:])
		for i in range(1, DL): ax.plot(x, result.iloc[i,:])

		# Shrinking the box size to append legend
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

		# Setting other plot properties
		ax.set(title="Trace of the nonlinear functions", xlabel="x", ylabel="F(x)")
		ax.grid(ls='--', c='lightgray')

		return True



	@staticmethod
	def apply_glasso_EM(**args):
		"""Implements EM Algorithm for group lasso

		Expected Args:
			beta_init : Initial Value for Beta
			A: Sufficient Statistics, a matrix 
			B: Sufficient Statistics
			ngroup: Number of Groups 
			lambda_value: LASSO Penalty 
			alpha2: EM decomposition parameter, tao^2 as mentioned in the paper 

		Returns:
			A Dictionary containing the following items:-
			sucess: Indictor for EM algorithm success
			beta: Updated Beta coefficients
		"""
		# Handling missing arguments
		required_args = ['beta_init', 'A', 'B', 'ngroup', 'lambda_value', 'alpha2']
		if set(required_args).intersection(set(args)) != set(required_args):
			print('[ERROR] Expected 6 required arguments in SLANTS.apply_glasso_EM(): \
				\n\t(1) beta_init = Initial Value for Beta \
				\n\t(2) A = Sufficient Statistics, a matrix \
				\n\t(3) B = Sufficient Statistics\
				\n\t(4) ngroup = Number of Groups \
				\n\t(5) lambda_value = LASSO Penalty \
				\n\t(6) alpha2 = EM decomposition parameter, tao^2 as mentioned in the paper')
			return False
		beta_init, A, B, ngroup, lambda_value, alpha2 = (args[element] for element in required_args)
		if beta_init.ndim != 2: beta_init = np.array([beta_init]) 

		# Convergence + Debug Parameters Recommended Setting
		K = 20 if 'K' not in args else args['K']
		tolerance_EM = 100 if 'tolerance_EM' not in args else args['tolerance_EM']
		eps_EM_abs = 1e-3 if 'eps_EM_abs' not in args else args['eps_EM_abs']
		eps_EM_rel = 1e-2 if 'eps_EM_rel' not in args else args['eps_EM_rel']
		debug_print = 0 if 'debug_print' not in args else args['debug_print'] # 0 => no print, 1 => show print, 2 => show print + beta plot

		p = np.size(beta_init)
		gsize = p//ngroup

		# Initialize the history container
		beta_history = np.copy(beta_init)
		beta = np.copy(beta_init)  

		######################## PHASE 1
		# Iterate between old beta -> compute r -> shrink r (new beta)
		for k in range(K):
				beta_old = np.copy(beta)
				r = beta.T - alpha2 * np.dot(A, beta.T) + np.dot(alpha2, B.T) # should be numpy.ndarray, need to check

				for j in range(ngroup):
						group = list(range((j*gsize), ((j+1) * gsize)))
						temp = alpha2 * lambda_value / np.sqrt(sum(np.square(r[group]))) # only works when r is np.ndarray

						if(temp < 1):
								beta[0, group] = (1-temp) * r[group].reshape(-1)
						else:
								beta[0, group] = 0
				beta_history = np.concatenate((beta_history, beta))

				# Checking convergence
				chg = np.sum(abs(beta - beta_old))
				if (chg < eps_EM_abs):
						if (debug_print):
								print("EM converges in {:d} iterations (abs)".format(k+1))
						break
				elif ((np.sum(abs(beta_old)) > 0) and (chg/np.sum(abs(beta_old)) < eps_EM_rel)): # /0 error rectified 
						if (debug_print):
								print("EM converges in {:d} iterations (rel)".format(k+1))
						break

		######################## PHASE 2
		# Post processing + checking convergence
		if ((np.sum(abs(beta_init))>0) & (np.sum(abs(beta - beta_init)) > tolerance_EM * np.sum(abs(beta_init)))):
				success = False
				print("EM blow up with rate {:.4f}".format(np.sum(abs(beta - beta_init))/np.sum(abs(beta_init))))
		else:
				success = True

		# Plot the history of first component
		if (debug_print > 1):
				length = np.size(beta_history[:,1], 0)
				x = np.arange(0, length)
				fig, ax = plt.subplots()
				ax.plot(x,beta_history[:, 1])
				fig.suptitle('EM history of beta[1]', fontsize=12)
				ax.set_xlabel('Iteration', fontsize=10) ; ax.set_ylabel("beta[1]",fontsize='medium')
				ax.xaxis.label.set_size(20)
				plt.draw()

		result = {
			'success':success, 
			'beta': beta
		}

		return result

	def main():
		print('==========< SLANTS Implementation >=============')
		print('This module packages the backbone structure for the Sequential Learning Algorithm for Nonlinear Time Series (SLANTS).')

	if __name__ == '__main__':
		main()