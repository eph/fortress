from fortress import SMCDriver

import os

os.system('export FC=gfortran')



from FRUIT import *
test_modules = ['test/test_random.f90',
               'test/test_prior.f90',
               'test/test_model.f90',
               'test/test_linalg.f90',
               'test/test_util.f90',
               'test/test_gensys.f90',
               'test/test_smc.f90',
               'test/test_particles.f90',
               'test/test_particle_filter.f90'
]

suite = test_suite(test_modules)
suite.build_run('test_driver.f90', 'make test_driver')
suite.summary()
