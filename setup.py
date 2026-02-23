#/urs/bin/env python

__version__ = "0.0.1-dev"

from distutils.core import setup

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Statitics
    Programming Language :: Python
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3.4
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""

setup(name='eurydice',
      version="0.0.1-dev",
      license='BSD3-AI',
      description="Python helper functions for ECHO microbiome data",
      long_description=("Library for handling a metadata data dictionary "
                        "(and breaking the fourth wall)."),
      author="J W Debelius",
      author_email="justine.debelius@jhu.edu",
      # maintainer="J W Debelius",
      # maintainer_email="jdebelius@ucsd.edu",
      packages=['eurydice',
		'eurydice.impute',
		'eurydice.impute.tests'
		'eurydice.plot'
	       ],
      install_requires=['numpy >= 1.24.4',
			'pandas >= 1.5.3',
			'seaborn >= 0.12.2',
			'matplotlib >= 3.6.0',
			'scipy >= 1.10.0',
			'scikit-bio >= 0.5.8'
                        'nose >= 1.3.7',
                        ],
      )
