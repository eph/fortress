import numpy as np
import os
import subprocess
import json
import tqdm
import pandas as p
import glob
my_env = os.environ.copy()
my_env['OPENBLAS_NUM_THREADS'] = '1'

if 'CONDA_PREFIX' in my_env.keys():
    inc_path = '$(CONDA_PREFIX)/include'
    lib_path = '$(CONDA_PREFIX)/lib'
elif 'ANACONDA_PATH' in my_env.keys():
    inc_path = '$(ANACONDA_PATH)/../include'
    lib_path = '$(ANACONDA_PATH)/../lib'
else:
    print('Error Setting PATHs for lib and including -- specify manually')
    inc_path = '$(ANACONDA_PATH)/../include'
    lib_path = '$(ANACONDA_PATH)/../lib'


def load_estimates(file_string, resample=True, paranames=None, posterior='final'):

    output_files = glob.glob(file_string)

    results = []
    for f in output_files:
        output_json = json.loads(open(f).read())

        posteriors = sorted([k for k in output_json.keys() if k.startswith('posterio')])

        if posterior=='final': 
            to_load = posteriors[-1]
        else:
            to_load = posterior

        res = p.DataFrame(output_json[to_load])
      
        if resample:
            inds = np.random.choice(res.shape[0], size=res.shape[0], p=res.weights)
            res = res.iloc[inds].reset_index()
      
        if paranames is not None:
            vs = [c for c in res.columns if c.startswith('var')]
            res = res.rename(columns=dict(zip(vs,paranames)))
          
        res['logmdd'] = np.array(output_json['Z_estimates']).sum()
      
        results.append(res)
    if len(output_files) == 1:
        return results[0]
    else:
        return p.concat(results, axis=0, keys=output_files)
    
    

#my_env['LD_LIBRARY_PATH']='/home/eherbst/Dropbox/code/fortress'

class SMCDriver(object):

    def __init__(self, executable):
        self.executable = executable
        res = subprocess.run([self.executable, '--help'],
                             env=my_env,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        self.help = res.stderr.decode(encoding='UTF-8')

    def run(self, **kwargs): 
        nproc = kwargs.pop('nproc',1)
        envs = kwargs.pop('env', {})
        output_file = kwargs.get('output_file','output.json')

        
        my_env.update(envs)

        mpi = 'mpirun -n {} '.format(nproc)
        args = sum([['--'+k.replace('_','-'),str(v)] for k,v in kwargs.items()], [])
        args = [a for a in args if a is not 'True']

        proc = subprocess.Popen([mpi+self.executable+' '+' '.join(args)],
                                env=my_env, stdout=subprocess.PIPE, shell=True,
                                stderr=subprocess.PIPE, universal_newlines=True) 

        pbar = tqdm.tqdm(total=1.0)
        i = 0
        for line in iter(proc.stdout.readline, ''):
            line2 = line.strip()
            if line2.startswith('iteration'):
                lab, it, of, tot = line2.split()
                pbar.update(1/float(tot))
                i = i + 1

        pbar.close()

        return json.loads(open(output_file).read())



makefile = """
LIB={lib_path}
INC={inc_path}
FPP=fypp
FRUIT=-I$(INC)/fruit -L$(LIB) -lfruit -Wl,-rpath=$(LIB)
FLAP=-I$(INC)/flap -L$(LIB) -lflap
FORTRESS=-I$(INC)/fortress -L$(LIB)/fortress -lfortress
JSON=-I$(INC)/json-fortran -L$(LIB)/json-fortran -ljsonfortran

FC={f90} -O3 -ffree-line-length-1000 #-Wall -fcheck=all -g -fbacktrace

smc_driver : {model_file} smc_driver.f90 
\t$(FC) $^  -I. -Wl,--start-group $(FORTRESS) $(JSON) $(FLAP) $(FRUIT) -l{lapack} -Wl,--end-group -o smc 

check_likelihood : {model_file} check_likelihood.f90 
\t$(FC) $^  -I. -Wl,--start-group $(FORTRESS) $(JSON) $(FLAP) $(FRUIT) -l{lapack} -Wl,--end-group -o check_likelihood 
"""

import os

driverfile ="""program smc_driver
  use iso_fortran_env, only: wp => real64

  use fortress, only: fortress_smc
  use model_t, only: model

  implicit none
  include 'mpif.h'

  type(fortress_smc) :: smc
  type(model) :: smc_model

  integer :: mpierror, rank, nproc, i

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

  smc_model = model()
  smc = fortress_smc(smc_model, nproc)
  call smc%estimate(rank)
  call mpi_finalize(mpierror)
end program smc_driver"""

checkfile ="""program check_likelihood
  use iso_fortran_env, only: wp => real64

  use model_t, only: model

  implicit none
  include 'mpif.h'

  type(model) :: smc_model

  real(wp) :: lik0 
  smc_model = model()
  lik0 = smc_model%lik(smc_model%p0)
  print*,'Likelihood @ p0: ', lik0
end program check_likelihood"""

simplefile ="""
module model_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: model_prior => prior 

  {other_includes}

  implicit none

  type, public, extends(fortress_abstract_bayesian_model) :: model

  contains
    procedure :: lik
  end type model

  interface model
    module procedure new_model
  end interface model

contains

type(model) function new_model() result(self)

character(len=144) :: name, datafile, priorfile
integer :: nobs, T, npara

name = 'qar'
datafile = 'data.txt'
priorfile = 'prior.txt'

nobs = 1
T = {T}
npara = {npara}

call self%construct_model(name, datafile, npara, nobs, T)

allocate(self%prior, source=model_prior(priorfile))

end function new_model

function lik(self, para, T) result(l)
class(model), intent(inout) :: self

real(wp), intent(in) :: para(self%npara)
integer, intent(in), optional :: T
real(wp) :: l

{lik}

end function lik

{other_functions}

end module model_t
"""

def make_model_file(lik,npara,T,other_functions='',other_includes=''):

    return simplefile.format(lik=lik, npara=npara, T=T, other_functions=other_functions,other_includes=other_inclues)


def make_smc(model_file, output_directory='_fortress_tmp', other_files=None,
             lib_path=lib_path,inc_path=inc_path, f90='mpif90', lapack='openblas',check=True):
    """
    Makes an smc driver.
    """
    tmp_dir = output_directory
    try:
        os.makedirs(tmp_dir)
    except FileExistsError:
        print('Directory exists!')
      
    with open(os.path.join(tmp_dir, 'model_t.f90'), 'w') as f:
        modelfile = model_file.format(output_directory=os.path.abspath(output_directory))
        for name, contents in other_files.items():
            basename = os.path.basename(name)
            outname = os.path.join(os.path.abspath(tmp_dir), basename)
            print('Replacing {} with {}'.format(basename, outname))
            modelfile = modelfile.replace(basename, outname)

        f.write(modelfile)

    with open(os.path.join(tmp_dir, 'makefile'), 'w') as f:
        f.write(makefile.format(model_file='model_t.f90',
                                lib_path=lib_path,inc_path=inc_path,f90=f90,lapack=lapack))

    with open(os.path.join(tmp_dir, 'smc_driver.f90'), 'w') as f:
        f.write(driverfile)

    with open(os.path.join(tmp_dir, 'check_likelihood.f90'), 'w') as f:
        f.write(checkfile)

    if other_files is not None:
      
        if isinstance(other_files, list):
            other_files = {f:f for f in other_files}

        for name, contents in other_files.items():

            basename = os.path.basename(name)

            outname = os.path.join(tmp_dir, basename)
            if isinstance(contents, (np.ndarray, p.DataFrame)):
                np.savetxt(outname, contents)
            elif os.path.isfile(contents):
                from shutil import copyfile
                copyfile(contents, outname)
            else:
                with open(outname, 'w') as f:
                    f.write(contents)



    proc = subprocess.run('cd {} && make smc_driver'.format(os.path.abspath(tmp_dir)),
                          shell=True, check=check, stderr=subprocess.STDOUT, env=my_env)

    return(SMCDriver(os.path.abspath(os.path.join(tmp_dir,'smc'))))
                
# def make_mcmc_driver(output_directory, model_file, other_files=None):
#     """
#     Makes an mcmc driver.
#     """
#     pass 
    
    
