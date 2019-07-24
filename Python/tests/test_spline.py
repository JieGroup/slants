"""Curates test cases by usng R's spline functions and tests Python code's results with the same"""

import pytest
import numpy as np
from numpy import testing as npt
from slants.model.spline import BSpline


# Test Cases for Checking if spline functions are working fine
# Note: order = k
tcase1_aug_knots = np.array([[0. ,0. ,0. ,0.125 ,0.25 ,0.375 ,0.5 ,0.625 ,0.75 ,0.875 , 1. , 1. , 1. ]])
@pytest.mark.parametrize(
	'boundary_sf, nBspline_sf, order_sf, augmented_knots_sf,\
	tt_sf, xx_sf, regressor_scalar_sf, \
	k_sf, x_sf, end_knot_sf, pos_poly_num_sf', # args in order of: augknt(), spval(), positive_polynomial()
	[
		([0,1], 10, 3, tcase1_aug_knots, # augknt()
		0.2, np.array([0, 0, 0, 0.5]), 0.36, #spval()
		3, 0, 1, 0.04) # positive_polynomial()
	]
)
def test_spline_func(boundary_sf, nBspline_sf, order_sf, augmented_knots_sf, tt_sf, xx_sf, regressor_scalar_sf, k_sf, x_sf, end_knot_sf, pos_poly_num_sf):
	'''Tests the validity of augknt(), spval() and positive_polynomial()
	Test Case Index: sf 

	Purpose: Check if the spline functions produce desired correct output
	'''
	npt.assert_array_almost_equal(BSpline.augknt(boundary=boundary_sf, nBspline=nBspline_sf, order=order_sf), augmented_knots_sf)
	assert pytest.approx(BSpline.spval(t=tt_sf, x=xx_sf)) == regressor_scalar_sf
	assert pytest.approx(BSpline.positive_polynomial(k=k_sf, nBspline=nBspline_sf, t=tt_sf, x=x_sf, end_knot=end_knot_sf)) == pos_poly_num_sf