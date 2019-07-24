import pytest
import numpy as np
from numpy import testing as npt
from slants.model.spline import BSpline


# Test Cases for augknt()
knotBox, nBspline, order = np.array([1.290209, 2.504524]), 10, 4 # knotbox is the max/min% quantile
@pytest.mark.parametrize(
	'args_augknt_fargs, result_augknt_fargs',
	[
		({'boundary': knotBox, 'nBspline': nBspline, 'order': order}, np.ndarray), # ==3 args
		({'boundary': knotBox, 'nBspline': nBspline, 'unreq_arg1': 1}, bool),
		({'boundary': knotBox, 'nBspline': nBspline, 'order': order, 'unreq_arg1': 1}, np.ndarray), # > 3 args
		({'boundary': knotBox, 'nBspline': nBspline, 'unreq_arg1': 1, 'unreq_arg2': 2}, bool),
		({'boundary': knotBox, 'nBspline': nBspline}, bool), # < 3 args
		({'boundary': knotBox, 'unreq_arg1': 1}, bool),
		({}, bool),
	]
)
def test_sufficient_function_args_augknt(args_augknt_fargs, result_augknt_fargs):
	'''Tests if augknt() can handle the variadic input arguments correctly
	Test Case Index: fargs 

	Purpose: Ensure future-proof compatibility of augknt() by accepting additional unrequired arguments

	Setup:
		# (req args) = Minimum Number of required arguments to be satisfied
		# (input args) = Number of input (variadic) arguments to the function 
	
	Range of Inputs:
		{Required Arguments}
			knotBox: Numpy array
			nBspline: Number of BSplines, +ve scalar
			order: Order of BSplines, +ve scalar
		{3 Unrequired Arguments: +ve scalars}

	Edge Cases: When # (req args) = 3;  
		(1) # (input args) < 3
		(2) # (input args) > 3
		(3) # (input args) == 3

	Positive Output: Numpy array Class type pertaining to array containing the knots, when # (req args) satisfied
	Negative Output: Bool Class type pertaining to False, when # (req args) not satisifed
	'''
	assert type(BSpline.augknt(**args_augknt_fargs)) == result_augknt_fargs



# Test Cases for positive_polynomial()
tt, x, k, end_knot, d = 0.2, 0, 3, 2, 0
@pytest.mark.parametrize(
	'args_poly_fargs, result_poly_fargs',
	[
		({'t': tt, 'x': x, 'k': k, 'end_knot': end_knot}, float), # ==4 args
		({'t': tt, 'x': x, 'k': k, 'd': 0}, bool),
		({'t': tt, 'x': x, 'k': k, 'end_knot': end_knot, 'd': 0}, np.float), # > 4 args
		({'t': tt, 'x': x, 'k': k, 'd': 0, 'unreq_arg1': 1}, bool),
		({'t': tt, 'x': x, 'k': k}, bool), # < 4 args
		({'t': tt, 'x': x, 'd': 0}, bool),
		({}, bool),
	]
)
def test_sufficient_function_args_poly(args_poly_fargs, result_poly_fargs):
	'''Tests if positive_polynomial() can handle the variadic input arguments correctly
	Test Case Index: fargs 

	Purpose: Ensure future-proof compatibility of positive_polynomial() by accepting additional unrequired arguments

	Setup:
		# (req args) = Minimum Number of required arguments to be satisfied
		# (input args) = Number of input (variadic) arguments to the function 
	
	Range of Inputs:
		{Required Arguments}
			t: Input value, scalar
			x: Knot value, scalar
			k: Order, +ve scalar
			end_knot: Largest value of knot, +ve scalar
		{Unrequired Arguments}
			d: 'dth' derivative, scalar

	Edge Cases: When # (req args) = 4;  
		(1) # (input args) < 4
		(2) # (input args) > 4
		(3) # (input args) == 4

	Positive Output = Flaat array Class type pertaining to +ve polynomial scalar, when # (req args) satisfied
	Negative Output = Bool Class type pertaining to False, when # (req args) not satisifed			
	'''
	
	assert type(BSpline.positive_polynomial(**args_poly_fargs)) == result_poly_fargs