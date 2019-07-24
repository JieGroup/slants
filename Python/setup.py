# Setup Information for the package
from setuptools import setup, find_packages
from os.path import dirname, join, realpath

with open(join(dirname(realpath(__file__)), "README.md")) as f:
	long_description = f.read()

setup(name='slants',
	version='1.0.0',
	description='Sequential Learning Algorithm for Nonlinear Time Series (SLANTS)',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/JieGroup/slants',
	author='Rachit Jas, Jiaqi Liu, Qiuyi Han, Jie Ding',
	author_email='dingj@umn.edu',
	license='GNU GPLv3',
	classifiers=[ 
		'Topic :: Education',
		'Topic :: Scientific/Engineering :: Information Analysis',
		'Topic :: Scientific/Engineering :: Mathematics',
		'Topic :: Scientific/Engineering',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',
		'Programming Language :: Python :: 3.8',
		'Intended Audience :: Education',
		'Intended Audience :: Developers',
		'Intended Audience :: Information Technology',
		'Intended Audience :: Science/Research',
	],
	keywords='sequential adaptive online time series machine learning expectation maximization algorithm group lasso statistics',
	packages=find_packages(exclude=['docs', 'tests*']),
	# package_dir={'slants': 'slants/'},
	package_data={'slants': ['data/*.csv']},
	install_requires=['numpy', 'pandas', 'matplotlib'],
)