"""SPLINE UTILITY FUNCTIONS
Low level B-spline approximation functions. All spline notations are consistent with those used here: http://www.stem2.org/je/bspline.pdf

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


# Importing Dependencies 
import numpy as np
from math import lgamma, exp, factorial

class BSpline:
	"""Utility class for configuring splines for SLANTS algorithm

	  	Description: 
	  		Defines B-spline approximation using power truncated polynomials. A closer look reveals the attempt to replicate MATLAB's spval and augknt spline functions  		
	"""

	@staticmethod
	def augknt(**args):
	  	"""Finds the knots for the B-Spline

	  	Expected Args:
	  		boundary: Boundary values of knots, defined as a list
	  		nBspline: Number of B-Splines
	  		order: Order of polynomial

	  	Returns:
	  		A numpy array containing the knots

	  	"""
	  	#print('augknt')		
  	
  		# Handling missing arguments
  		required_args = ['boundary', 'nBspline', 'order']
	  	if set(required_args).intersection(set(args)) != set(required_args):
	  		print('[ERROR] Expected 3 required arguments in BSpline.augknt(): \
	  			\n\t(1) boundary = Boundary of knots \
	  			\n\t(2) nBspline = Number of B-Splines\
	  			\n\t(3) order = Order of polynomoial')
	  		return False

	  	# Returning the knots for the B-Spline
	  	return np.hstack((np.repeat(args['boundary'][0], args['order']-1), \
	  		np.linspace(args['boundary'][0], args['boundary'][1], args['nBspline']-args['order']+2), \
	  		np.repeat(args['boundary'][1], args['order']-1))).reshape(1, -1)


	@staticmethod
	def positive_polynomial(**args):
	  	"""Returns a positive-polynomial scalar

	  	Expected Args:
	  		t: Input Value, scalar value
	  		x: Knot Value, scalar value
	  		k: Order
	  		end_knot: Largest value for knot
	  		d: 'dth' derivative, if required

	  	Returns:
	  		A number of positive polynomial result

	  	Example: 
	  		positive_polynomial(t=0.2, x=0, k=3, end_knot=2, d=0)
	  	"""		
		#print('positive_polynomial')

  		# Handling missing arguments
  		required_args = ['t', 'x', 'k', 'end_knot']
	  	if set(required_args).intersection(set(args)) != set(required_args):	  		
	  		print('[ERROR] Expected 4 required arguments in BSpline.positive_polynomial(): \
	  			\n\t(1) t = Input Value \
	  			\n\t(2) x = Knots Value\
	  			\n\t(3) k = Order\
	  			\n\t(4) end_knot = Largest value for knot')
	  		return False
	  	d = 0 if 'd' not in args else args['d']

	  	# Storing required arguments
	  	t, x, k, end_knot = args['t'], args['x'], args['k'], args['end_knot']

	  	# Returning positive polynomial scalar if input is correct + d > 0
	  	if (t <= x) or (d > k-1) or (t > end_knot): return 0
	  	elif (d > 0): return exp(lgamma(k) - lgamma(k-d)) * (t - x)**(k - 1 - d)
	  	else: return (t - x)**(k - 1)


	@staticmethod
	def spval(**args):  
	  	"""Calculates the divided difference and thus a single value of regressor

	  	Description: 
	  		Takes the divided difference with order k = length(x)-1
	  		x should be atleast length 2, and in non-decreasing order
	  		Duplicate values in x are allowed 

	  	Expected Args:
	  		t: Input Value, and is a scalar
	  		x: Knots, and is a vector of length 'k + 1'

	  	Returns:
	  		Regressor scalar value
	  	"""		
		#print('spval')

  		# Handling missing arguments
	  	required_args = ['t', 'x']
	  	if set(required_args).intersection(set(args)) != set(required_args):
	  		print('[ERROR] Expected 2 required arguments in BSpline.spval(): \
	  			\n\t(1) t = Input Value, and is a scalar \
	  			\n\t(2) x = Knots, and is a vector of length \'k + 1\'')
	  		return False

	  	# Storing required arguments
	  	t, x = args['t'], args['x']

	  	# Defining the order and end of the knot
	  	k, end_knot = len(x) - 1, max(x)
	  	phi = np.apply_along_axis(lambda x: BSpline.positive_polynomial(t=t, x=x, k=k, end_knot=end_knot), \
	  		axis=0, arr=x.reshape(1,-1)).reshape(-1)

	  	# Computing 'phi'
	  	for i in range(1,k+1):
	  		phi_new = np.repeat(np.nan, k+1-i)

	  		for j in range(k+1-i):
	  			if x[j] == x[j+i]:
	  				phi_new[j] = abs(BSpline.positive_polynomial(t=t, x=x[j], k=k, end_knot=end_knot, d=i))/factorial(i)
	  			else:
	  				phi_new[j] = abs((phi[j+1] - phi[j]) / (x[j+i] - x[j]))

	  		phi = phi_new

	  	return np.asscalar((end_knot-x[0]) * phi)

	def main():
	  	print('==========< BSpline Implementation >=============')
	  	print('This module packages all spline utility and approximation functions.')

	if __name__ == '__main__':
	  	main()	  	