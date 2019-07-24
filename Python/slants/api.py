"""SLANTS API WRAPPER
Friendly user interface for accessing the SLANTS algorithm functions    

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
# Dummy code, TODO: Add Online Model
#import slants.model.base as bs
from slants.model.base import SLANTS
from os.path import exists, join, realpath
from os import listdir
from re import findall
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np


class Model:
	"""Utility Class for the Raw Implementation of the SLANTS algorithm

	Description:
		Provisions the underlying logic behind the building blocks of the algorithm, such as the data preprocessing and Glasso EM steps  
	""" 
	def fit(self, data, predict, maximum_lag, spline_info, forget_factor, lambda_init, tau2_init, sparse_lambda_tolerance, delta_multiplier, shrink_step_size, 
		online=True, debug_print=False, test_size=50, idle_tau2_tolerance=50, safe_shrink_lambda=10**(0.4), zero_coef_tolerance=3, shrink_tau2=1.1, random_beta=False):
		"""Fits the SLANTS model on data
			Note: Data is automatically prepared (standardized + transformed as in the D̃ = Dl step in the paper). Splines are also inherently handled

		Expected Args: {format - [type, description]}
			data: Numpy array; Original data matrix, X
			predict: Scalar; Column to indicate the random variable to be predicted
			maximum_lag: Scalar; Maximum lag, L, to be considered
			spline_info: Dictionary; initial spline configuration containing 'order' and 'num_Bsplines'
			forget_factor: Numpy vector; comprising ideally decreasing weights in range [0,1] indicating the importance of past values (generally, less significance for older values)
			lambda_init: Scalar; LASSO penalty initial value
			tau2_init: Scalar; EM decomposition parameter initial value
			sparse_lambda_tolerance: Scalar; Favors smaller lambda within spaTol_c tolerance
			delta_multiplier: Scalar; Small step size multiplier, to move along the lambda channels
			shrink_step_size: Numpy vector; Adjusts lambda before they get out of control, recommended to be equal to 'forget_factor' vector values in case of confusion 


		Optional Args: {format - [default value, description]}
			online: True, Boolean to indicate whether online version of SLANTS should be utilized using previous model parameters. Automatically handles situations wherein no past serialized model parameters exist 		
			debug_print: False, Boolean to indicate internal status of algorithm involving debug print statements
			test_size: 50, Number of observations to average prediction error
			idle_tau2_tolerance: 50, Increase tau^2 if it idles for a while
			safe_shrink_lambda: 10^(0.4), Adjust lambda before they get out of control
			zero_coef_tolerance: 3, Shrink lambda if coefficients are 0 for a long time
			shrink_tau2: 1.1, Shrink tau^2     
			random_beta: False, Use custom provided beta as initializers for the random coefficients

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
		try:

			arr, scalar, bools = [data, forget_factor, shrink_step_size], [predict, maximum_lag, lambda_init, tau2_init, sparse_lambda_tolerance, delta_multiplier, 
			test_size, idle_tau2_tolerance, safe_shrink_lambda, zero_coef_tolerance, shrink_tau2], [online, debug_print]

			for entry in arr: 
				if not isinstance(entry, (np.ndarray)): raise ValueError('Numpy array data type not satisfied.')
			for entry in scalar: 
				if not type(entry) in (float, int): raise ValueError('Scalar data type not satisfied.')
			for entry in bools: 
				if not isinstance(entry, bool): raise ValueError('Boolean data type not satisfied.')
			if not isinstance(spline_info, dict): raise ValueError('Dictionary data type not satisfied for spline configuration.')
			if not isinstance(random_beta, (np.ndarray, bool)): raise ValueError('Data type not satisfied for random beta coefficienits.')

			# Preprocess data
			preproc_results = SLANTS.preprocess_data(X=data, L=maximum_lag, order=spline_info['order'], nBspline=spline_info['num_Bsplines'])
			transformed_data, spline_config = preproc_results['x'], preproc_results['spconfig']
			print('==========< Data Preparation Finished >=============')

			fit_results = SLANTS.fit(y=data[:, predict-1], x=transformed_data, 
				D=data.shape[1], L=maximum_lag, lambda_value=forget_factor, gamma_init=lambda_init, alpha2_init=tau2_init, spconfig=spline_config, 
				move_size=delta_multiplier, spaTol_gamma=sparse_lambda_tolerance, shrink_step_size=shrink_step_size, y_index=predict-1, test_size=test_size, 
				tol_idle_alpha2=idle_tau2_tolerance, safe_shrink_gamma=safe_shrink_lambda, tol_all0times=zero_coef_tolerance, shrink_alpha2=shrink_tau2, 
				online=online, debug_print=debug_print, use_random_beta=random_beta)

			self._spline_config, self._model_params = spline_config, fit_results
		except ValueError as e:
			print('[ERROR] A <ValueError> exception occured in Model.fit(). \nExplanation: {} \nTip: Check data type or relevance of model parameter values, like valid maximum lag values and spline info configurations.'.format(e.args[0]))						
			return False			
		except Exception as e:
			print('[ERROR] Oops! An unexpected error occured, of type: {} in Model.fit()'.format(type(e)))
			return False			

	def model_params(self, filename='None'):
		"""Returns a dictionary containing the model parameters 

		Optiona; Args:
			filename: Name of the serialized model file 

		Returns:
			A Dictionary containing the following items:-
			prediction_err: Historical predicted error; 1st column indicates 1st channel, etc. 
			prediction_val: Historical predicted values; 1st column indicates 1st channel, etc. 
			lambda_history: Historical lambda channels; 1st channel is the smallest one, 3rd channel is the largest one   
			lambda_opt: Optimized lambda; 2nd channel in lambda_history 
			tau2_opt: Optimized Tau^2 
			coef_opt: Best coefficient, the last row and second channel related in coef_history  
			coef_history: All stored coefficient
			spline_config: Spline configuration
			sufficient_stats_A, sufficient_stats_B: Sufficient statistics
		"""
		# Checking for model parameter existence errors
		try:
			rename_list = [('gamma_history', 'lambda_history'), ('gamma_opt', 'lambda_opt'), ('alpha_opt', 'tau2_opt'),
							('spconfig', 'spline_config'), ('pre_err', 'prediction_err'), ('predict1', 'prediction_val'),
							('beta_history', 'coef_history'), ('beta_opt', 'coef_opt'), ('mY', 'mean_of_predict_rand_var'),
							('A', 'sufficient_stats_A'), ('B', 'sufficient_stats_B')]
			if  filename != 'None': 
				if exists(filename):
					with open(filename, 'rb') as serde_handle:
						model_params = pickle.load(serde_handle)
						for name in rename_list: model_params[name[1]] = model_params.pop(name[0]) 
						return model_params
				else:
					print('[ERROR] Serialized model file not found. Tip: Check for filename existence or spelling errors.')
					return False

			if not hasattr(self, '_model_params'):
				print('[ERROR] The SLANTS model should be fit first. Tip: Use fit() method.')
				return False

			model_params = deepcopy(self._model_params)
			for name in rename_list: model_params[name[1]] = model_params.pop(name[0])
			return model_params
		except Exception as e:
			print('[ERROR] Oops! An unexpected error occured, of type: {} in Model.model_params()'.format(type(e)))
			return False


	def plot(self, kind='coef'):
		"""Graphs a type of plot associated with the model to track the tuning parameters  

		Expected Args:
			kind: Type of requested plot to be graphed using the model parameters. At the moment, 4 types of plot can be created:
				(1) lambda = Optimal Lambda Plot
				(2) tau = Optimal Tau Plot
				(3) coef = Coefficient-Time Plot
				(4) trace = Trace of the nonlinear functions


		Returns:
			False, if a type mismatch occurs		
		"""
		# Checking for model parameters existence
		try:
			if not hasattr(self, '_model_params'):
				print('[ERROR] The SLANTS model should be fit first. Tip: Use fit() method.')
				return False

			if kind == 'trace':
				status = Plotter.plot_tuning_params(opt_params=self._model_params['beta_opt'], supp_params=self._spline_config, kind=kind)
			elif kind == 'lambda':
				status = Plotter.plot_tuning_params(opt_params=self._model_params['gamma_opt'], supp_params='N/A', kind=kind)
			elif kind == 'tau':
				status = Plotter.plot_tuning_params(opt_params=self._model_params['alpha_opt'], supp_params='N/A', kind=kind)  			
			else:
				status = Plotter.plot_tuning_params(opt_params=self._model_params['beta_history'], supp_params='N/A', kind=kind)			

			return True if status else False
		except Exception as e:
			print('[ERROR] Oops! An unexpected error occured, of type: {} in Model.plot()'.format(type(e)))
			return False

	def reset_model(self):
		"""Reset internal model parameters and spline configuration, maybe useful if the model is not performing as required   

		Expected Args:
			None

		Returns:
			False, if SLANTS model is already in reset mode 		
		"""
		# Checking for model parameters existence
		try:
			if not hasattr(self, '_model_params') or not hasattr(self, '_spline_config'):
				print('[ERROR] The SLANTS model is not required to be reset. Tip: Use fit() method to compute model parameters.')
				return False

			del self._model_params
			del self._spline_config
			return True
		except Exception as e:
			print('[ERROR] Oops! An expected error occured, of type: {} in Model.reset_model()'.format(type(e)))
			return False


	def help(self, method='fit', ls=False):
		"""Wrapper to show a Linux 'man'-like help documentation for methods/class   

		Expected Args:
			method: Method name whose documentation is required. Possible values:
				(1) fit = Fit SLANTS model
				(2) model_params = Retrieve SLANTS model parameters
				(3) plot = View diagnostic plots
				(4) reset_model = Reset model in case of problems
				(5) help = Display information of help method
				(6) all = Display entire documentation of class 
			ls: Shows the list of available methods 
		"""
		try:
			if ls:
				print('Here are the list of available methods: \n\t {}'.format(' | '.join((attribs for attribs in dir(self) if '__' not in attribs and 
					'_model_params' not in attribs and '_spline_config' not in attribs))))
				return True

			if method not in set(['fit', 'model_params', 'plot', 'reset_model', 'help', 'all']):
				print('<WARNING> Looks like "{}"" method does not exist. Here are the list of available methods: \n\t {}'.format(method, ' | '.join((attribs for attribs in dir(self) if '__' not in attribs and 
					'_model_params' not in attribs and '_spline_config' not in attribs))))
				if input('Were you looking to check the entire class documentation instead? (y/n)') == 'y':
					print('\n\n')
					help(self)
					return True
				else:
					print('<WARNING> Try again with one of the valid methods listed above.')
					return False
			
			if method != 'all':
				exec('help(self.{})'.format(method))
				return True		
			help(self)
			return True
		except Exception as e:
			print('[ERROR] Oops! An unexpected error occured, of type: {} in Model.help()'.format(type(e)))
			return False								


class Plotter:
	"""Utility Class for handling plots associated with the SLANTS algorithm

	Description:
		Encapsulates utility methods for the purpose of generating plots    
	""" 	 	
	@staticmethod
	def plot_tuning_params(**args):
		"""Handles graphing different requested plots

		Expected Args:
			opt_params: Optimal Model Parameters
			supp_params: Supporting Parameters for creating the plots
			kind: Type of requested plot to be graphed using the model parameters. At the moment, 4 types of plot can be created:
				(1) lambda = Optimal Lambda Plot
				(2) tau = Optimal Tau Plot
				(3) coef = Coefficient-Time Plot
				(4) trace = Trace of the nonlinear functions


		Returns:
			False, if a type mismatch happens
		"""  	
		# Handling missing arguments
		required_args = ['opt_params', 'supp_params', 'kind']
		if set(required_args).intersection(set(args)) != set(required_args):  		
			print('[ERROR] Expected 3 required arguments in Plotter.plot_tuning_params(): \
				\n\t(1) opt_params = Optimal Model Parameters \
				\n\t(2) supp_params = Supporting Parameters for creating the plots\
				\n\t(3) kind = Type of requested plot to be graphed using the model parameters')
			return False
		opt_params, supp_params, kind = (args[element] for element in required_args)
		if kind.lower() not in set(['lambda', 'tau', 'coef', 'trace']):
			print('[ERROR] Expected one of these 4 plot types in Plotter.plot_tuning_params(): \
				\n\t(1) lambda = Optimal Lambda Plot \
				\n\t(2) tau = Optimal Tau Plot \
				\n\t(3) coef = Coefficient-Time Plot \
				\n\t(4) trace = Trace of the nonlinear functions')
			return False

		if kind == "trace":
			return SLANTS.plot_coeff(beta=opt_params, knots=supp_params['knots'], nBspline=supp_params['nBspline'], order=supp_params['order'])

		# Setting up figure and subplot
		fig = plt.figure()
		ax = plt.subplot(111)

		# Plotting figures
		if kind in ('lambda', 'tau'):
			opt_params_symb = 'λ' if kind == 'lambda' else 'τ'  
			ax.plot(opt_params)
			ax.set(title="Optimal {}".format(opt_params_symb), xlabel="t", ylabel=opt_params_symb)
		else:
			mid_channel_index = opt_params.shape[1]//3			
			for i in range(0, mid_channel_index):
				ax.plot(opt_params[:,mid_channel_index+i])
				ax.set(title="All {} coefficients with time".format(mid_channel_index), xlabel="Time", ylabel="Coefficient Values")
		ax.grid(ls='--', c='lightgray')

		return True

class Experiment:
	"""Utility Class for providing sample experiments for testing/tutorial purposes

	Description:
		Easy-to-use API for accessing sample experiment data    
	"""
	@staticmethod
	def load(experiment_num):
		"""Loads data from CSV synthetic experiment file

		Expected Args:
			experiment_num: Number ID of the synthetic experiment

		Returns:
			False, if a type mismatch happens
		""" 
		try:
			inbuilt_exps = 3
			if experiment_num == '?':
				print('{} Synthetic experiments available at this time. Tip: Try using a number in the range [1, {}] as an argument to the Experiment.load() method.'.format(inbuilt_exps, inbuilt_exps))
				return True				
			if experiment_num not in list(range(1, inbuilt_exps+1)): #as further experiments are added
				print('[ERROR] Synthetic experiment {} not available at this time.'.format(experiment_num))
				return False

			# Initializing paths
			synth_exp_path_presuffix = ('synthexp{}_'.format(experiment_num), '.csv')
			base_dir = join(realpath(__file__)[:-6], "data")
			filenames = listdir(base_dir)
			fqpn_path = [join(base_dir, bdr) for bdr in listdir(base_dir)]
			# print(synth_exp_path_presuffix)
			# for x,y in zip(filenames, fqpn_path):
			# 	print(x, y)

			# Loading synthetic experiment data
			synth_exp_data = {}
			for file, path in zip(filenames, fqpn_path):
				# print(path, file, flush=True)
				valid = findall(r'{}(\w+){}'.format(synth_exp_path_presuffix[0], synth_exp_path_presuffix[1]), file)
				if valid: key = valid[0]
				else: continue
				# print(key)
				synth_exp_data[key] = np.genfromtxt(path, delimiter=',')
				
				# Handle single value parameters
				if synth_exp_data[key].size == 1:
					single_val = synth_exp_data[key].item()
					if isinstance(single_val, bool):
						synth_exp_data[key] = bool(single_val)
					if int(single_val) == float(single_val):
						synth_exp_data[key] = int(single_val)
					else:
						synth_exp_data[key] = float(single_val)				

			return synth_exp_data
		except Exception as e:
			print('[ERROR] Oops! An unexpected error occured, of type: {} in Experiments.load()'.format(type(e)))
			return False

		# print(getcwd(), flush=True)
		# print(abspath(__file__))
		# print('{}data'.format(realpath(__file__)[:-6]))
		# listdir('./data/')
		#@the.desert.eagle	