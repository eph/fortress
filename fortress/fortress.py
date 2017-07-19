import numpy as np
import os
import subprocess
import json
import tqdm
import pandas as p
my_env = os.environ.copy()
my_env['OPENBLAS_NUM_THREADS'] = '1'
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

        return json.loads(open('output.json').read())



makefile = """
LIB={lib_path}
INC={inc_path}
FPP=fypp
FRUIT=-I$(INC)/fruit -L$(LIB) -lfruit -Wl,-rpath=$(LIB)
FLAP=-I$(INC)/flap -L$(LIB) -lflap
FORTRESS=-I$(INC)/fortress -L$(LIB)/fortress -lfortress
JSON=-I$(INC)/json-fortran -L$(LIB)/json-fortran -ljsonfortran

FC={f90} -O3 #-Wall -fcheck=all -g -fbacktrace

smc_driver : smc_driver.f90 {model_file}
\t$(FC) $^  -I. -Wl,--start-group $(FORTRESS) $(JSON) $(FLAP) $(FRUIT) -l{lapack} -Wl,--end-group -o smc 
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

def make_smc(model_file, output_directory='_fortress_tmp', other_files=None,
             lib_path='$(HOME)/anaconda3/lib',inc_path='$(HOME)/anaconda3/include',
             f90='mpif90', lapack='openblas'):
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
        f.write(makefile.format(model_file='model_t.f90',lib_path=lib_path,inc_path=inc_path,f90=f90,lapack=lapack))

    with open(os.path.join(tmp_dir, 'smc_driver.f90'), 'w') as f:
        f.write(driverfile)

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
                          shell=True, check=True, stderr=subprocess.STDOUT, env=my_env)

    return(SMCDriver(os.path.abspath(os.path.join(tmp_dir,'smc'))))
                
# def make_mcmc_driver(output_directory, model_file, other_files=None):
#     """
#     Makes an mcmc driver.
#     """
#     pass 
    
    
