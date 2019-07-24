"""Curates test cases by usng R's model functions and tests Python code's results with the same"""

# Importing Dependencies
import pytest
import numpy as np
from numpy import testing as npt
from slants.model.base import SLANTS

# Helper function for loading test case results for preprocess_data test function
def tcase_load_res_pd(tcase_num):
	tcase_path_presuffix = ('./test_cases/base/get_regressor/tcase{}_R_'.format(tcase_num), '.csv') 
	fb = np.genfromtxt('{}feasibleBox{}'.format(*tcase_path_presuffix), delimiter=',')
	kb = np.genfromtxt('{}knotBox{}'.format(*tcase_path_presuffix), delimiter=',')
	k = np.genfromtxt('{}knots{}'.format(*tcase_path_presuffix), delimiter=',')
	sc = np.genfromtxt('{}scale{}'.format(*tcase_path_presuffix), delimiter=',')
	spkb = np.genfromtxt('{}spconfig_knotBox{}'.format(*tcase_path_presuffix), delimiter=',')
	spk = np.genfromtxt('{}spconfig_knots{}'.format(*tcase_path_presuffix), delimiter=',')
	spnb = np.genfromtxt('{}spconfig_nBspline{}'.format(*tcase_path_presuffix), delimiter=',')
	spo = np.genfromtxt('{}spconfig_order{}'.format(*tcase_path_presuffix), delimiter=',')
	x = np.genfromtxt('{}x{}'.format(*tcase_path_presuffix), delimiter=',')
	return dict(fb=fb, kb=kb, k=k, sc=sc, sp=dict(knotBox=spkb, knots=spk, nBspline=spnb, order=spo), x=x)


# Helper function for loading test case results for fit test function
fit_data_str_unprocessed_list = 'alpha_opt, alpha2_init, beta_history, beta_opt, D, feasibleBox, gamma_history, gamma_init, gamma_opt, \
knotBox, knots, L, lambda, move_size, mY, N, nBspline, order, pre_err, predict1, random_beta1, safe_shrink_gamma, \
scale, shrink_step_size, spaTol_gamma, spconfig_knotBox, spconfig_knots, spconfig_nBspline, spconfig_order, \
test_size, X, xx, y'.split(', ')

fit_data_list = [fit_data.strip() for fit_data in fit_data_str_unprocessed_list]

def tcase_load_data_fit(tcase_num):
	# TODO: compress the loader code if possible  
	tcase_path_presuffix = ('./test_cases/base/fit/tcase{}_fit_'.format(tcase_num), '.csv')
	fit_loaded_data = {}
	for fit_data in fit_data_list:
		fit_loaded_data[fit_data] = np.genfromtxt('{}{}{}'.format(tcase_path_presuffix[0], fit_data, tcase_path_presuffix[1]), delimiter=',')
		if fit_loaded_data[fit_data].size == 1:
			if int(fit_loaded_data[fit_data].item()) == float(fit_loaded_data[fit_data].item()):
				fit_loaded_data[fit_data] = int(fit_loaded_data[fit_data].item())
			else:
				fit_loaded_data[fit_data] = float(fit_loaded_data[fit_data].item())				
	return fit_loaded_data

# Test Cases for validating get_regressor() results
# Note: Index indicator is always an array of same column index to indicate which knot column is being used for computing regressors of that input vector 
tcase1_knots=np.array([[ 1.29020929,  4.45248695,  7.61476461, 10.77704227],
       [ 1.29020929,  4.45248695,  7.61476461, 10.77704227],
       [ 1.29020929,  4.45248695,  7.61476461, 10.77704227],
       [ 1.29020929,  4.45248695,  7.61476461, 10.77704227],
       [ 1.4636828 ,  4.62596046,  7.78823812, 10.95051578],
       [ 1.63715632,  4.79943398,  7.96171164, 11.1239893 ],
       [ 1.81062984,  4.9729075 ,  8.13518516, 11.29746282],
       [ 1.98410335,  5.14638101,  8.30865868, 11.47093634],
       [ 2.15757687,  5.31985453,  8.48213219, 11.64440985],
       [ 2.33105039,  5.49332805,  8.65560571, 11.81788337],
       [ 2.50452391,  5.66680157,  8.82907923, 11.99135689],
       [ 2.50452391,  5.66680157,  8.82907923, 11.99135689],
       [ 2.50452391,  5.66680157,  8.82907923, 11.99135689],
       [ 2.50452391,  5.66680157,  8.82907923, 11.99135689]])
tcase1_reg_vec=np.array([0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.1249184557,
0.6586540297, 0.2162992620, 0.0001282526, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000])
@pytest.mark.parametrize(
	'input_data_vec_gr, index_indicator_gr, spconfig_gr, regressors_vec_gr',
	[
		(
			np.linspace(1, 5, 5), np.repeat(0, 5), dict(order=4, nBspline=10, knots=tcase1_knots), tcase1_reg_vec
		)
	]
)
def test_get_regressor(input_data_vec_gr, index_indicator_gr, spconfig_gr, regressors_vec_gr):
	'''Tests the validity of get_regressor()
	Test Case Index: gr

	Purpose: Check if get_regressor() produces the desired correct output
	'''
	#print(SLANTS.get_regressor(x=input_data_vec_gr, d=index_indicator_gr, spconfig=spconfig_gr))
	npt.assert_array_almost_equal(SLANTS.get_regressor(x=input_data_vec_gr, d=index_indicator_gr, spconfig=spconfig_gr), regressors_vec_gr)



# Test Cases for validating preprocess_data() results
@pytest.mark.parametrize(
	'X_pd, L_pd, order_pd, nBspline_pd, results_pd',
	[
		(
			np.linspace(1,100,100).reshape(10,10,order='F'), 5, 4, 15, tcase_load_res_pd(1) # For X=1:100
		)
	]
)
def test_preprocess_data(X_pd, L_pd, order_pd, nBspline_pd, results_pd):
	'''Tests the validity of preprocess_data()
	Test Case Index: pd

	Purpose: Check if preprocess_data() produces the desired correct output
	'''	
	results = SLANTS.preprocess_data(X=X_pd, L=L_pd, order=order_pd, nBspline=nBspline_pd)
	npt.assert_array_almost_equal(results['feasibleBox'], results_pd['fb'])
	npt.assert_array_almost_equal(results['knotBox'], results_pd['kb'])
	npt.assert_array_almost_equal(results['knots'], results_pd['k'])
	npt.assert_array_almost_equal(results['scale'], results_pd['sc'])
	npt.assert_array_almost_equal(results['spconfig']['knotBox'], results_pd['sp']['knotBox'])
	npt.assert_array_almost_equal(results['spconfig']['knots'], results_pd['sp']['knots'])
	assert results['spconfig']['nBspline'] == results_pd['sp']['nBspline']
	assert results['spconfig']['order'] == results_pd['sp']['order']
	npt.assert_array_almost_equal(results['x'], results_pd['x'])

# Test Cases for validating apply_glasso_EM() results
@pytest.mark.parametrize(
	'beta_init_age, A_age, B_age, ngroup_age, lambda_value_age, alpha2_age, results_age', # beta, A and B => structured as 2D array each
	[
		(
			# np.matrix('1.0,2.0,3.0'), np.matrix('1.0,2.0,3.0;1.0,2.0,3.0;1.0,2.0,3.0'), np.matrix('1.0,2.0,3.0'), 3, 0.05, 0.05, dict(success=True, beta=[[-2.217133, -0.252133,  1.687867]]) 
			np.array([np.linspace(1,3,3)]), np.array([[1,2,3], [1,2,3], [1,2,3]]), np.array([np.linspace(1,3,3)]), 3, 0.05, 0.05, dict(success=True, beta=[[-2.217133, -0.252133,  1.687867]]) 
		)
	]
)
def test_apply_glasso_EM(beta_init_age, A_age, B_age, ngroup_age, lambda_value_age, alpha2_age, results_age):
	'''Tests the validity of test_apply_glasso_EM()
	Test Case Index: age

	Purpose: Check if apply_glasso_EM() produces the desired correct output
	'''	
	results = SLANTS.apply_glasso_EM(beta_init=beta_init_age, A=A_age, B=B_age, ngroup=ngroup_age, lambda_value=lambda_value_age, alpha2=alpha2_age)
	assert results['success'] == results_age['success']
	npt.assert_array_almost_equal(results['beta'], results_age['beta'])




# Deprecated test until further notice 
# Test Cases for validating fit() results
# fit_loaded_data1 = tcase_load_data_fit(1)
# tc1 = tuple([fit_loaded_data1[x] for x in fit_data_list])
# print(tc1)
#print(set(tuple([type(fit_loaded_data[x]) for x in fit_data_list])))
#print(len(tuple([x for x in fit_data_list])))


# @pytest.mark.parametrize(
# 	'alpha_opt, alpha2_init, beta_history, beta_opt, D, feasibleBox, gamma_history, gamma_init, gamma_opt, \
# 	knotBox, knots, L, lambda_value, move_size, mY, N, nBspline, order, pre_err, predict1, random_beta, safe_shrink_gamma, \
# 	scale, shrink_step_size, spaTol_gamma, spconfig_knotBox, spconfig_knots, spconfig_nBspline, spconfig_order, \
# 	test_size, X, xx, y',
# 	[
# 		tuple(tcase_load_data_fit(1)[x] for x in fit_data_list) #test case number: 1
# 	]
# )
# def test_fit(alpha_opt, alpha2_init, beta_history, beta_opt, D, feasibleBox, gamma_history, gamma_init, gamma_opt,
# 	knotBox, knots, L, lambda_value, move_size, mY, N, nBspline, order, pre_err, predict1, random_beta, safe_shrink_gamma,
# 	scale, shrink_step_size, spaTol_gamma, spconfig_knotBox, spconfig_knots, spconfig_nBspline, spconfig_order,
# 	test_size, X, xx, y):
	'''Tests the validity of test_apply_glasso_EM()
	Test Case Index: age

	Purpose: Check if apply_glasso_EM() produces the desired correct output
	'''	
	# spconfig = {'order': spconfig_order, 'nBspline': spconfig_nBspline, 'knots': spconfig_knots, 'knotBox': spconfig_knotBox}
	# results = SLANTS.preprocess_data(X=X, L=L, order=order, nBspline=nBspline)
	# npt.assert_array_almost_equal(results['feasibleBox'], feasibleBox)
	# npt.assert_array_almost_equal(results['knotBox'], knotBox)
	# npt.assert_array_almost_equal(results['knots'], knots)
	# npt.assert_array_almost_equal(results['scale'], scale)
	# npt.assert_array_almost_equal(results['spconfig']['knotBox'], spconfig['knotBox'])
	# npt.assert_array_almost_equal(results['spconfig']['knots'], spconfig['knots'])
	# assert results['spconfig']['nBspline'] == spconfig['nBspline']
	# assert results['spconfig']['order'] == spconfig['order']
	# npt.assert_array_almost_equal(results['x'], xx)
	#print(gamma_init)	
	# results = SLANTS.fit(
	# 	y=y, x=xx, D=D, L=L, lambda_value=lambda_value, gamma_init=gamma_init, alpha2_init=alpha2_init, spconfig=spconfig,
	# 	move_size=move_size, spaTol_gamma=spaTol_gamma, shrink_step_size=shrink_step_size, test_size=test_size, 
	# 	safe_shrink_gamma=safe_shrink_gamma, debug_print=1 
	# 					)