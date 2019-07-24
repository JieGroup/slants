# import pytest
# import numpy as np
# from numpy import testing as npt
# from slants.model.spline import BSpline
# from slants.model.base import SLANTS

# def test_func():
# 	# Set 1: Strict Float Comparison
# 	x = np.array([1., 1e-10, 1e-20])
# 	eps = np.finfo(x.dtype).eps
# 	print(x*eps/2 + x)	
# 	npt.assert_array_almost_equal_nulp(x, x*eps/2 + x)

# 	# Set 2: Solution to Strict Float Comparison
# 	print(x*eps + x)
# 	#npt.assert_array_almost_equal_nulp(x, x*eps + x) # this fails, more stronger comparison
# 	npt.assert_array_almost_equal(x, x*eps + x) # this won't, chiller comparison

# 	# Set 3: Simple np.array example  
# 	x=np.array([1,2,3,74.1])
# 	y=np.array([1,2,3,74.1000001499])	
# 	npt.assert_array_almost_equal(x,y,decimal=7)
# 	print(x, y)

# 	# Set 4: Checking if the package's module's classes are loaded correctly  
# 	print(BSpline.spval(t=1,x=np.array([2,3])))
# 	print(SLANTS.get_regressor())


