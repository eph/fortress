SRC=src
TEST=test
VPATH=.:$(SRC):$(TEST):templates

LIBOBJS=fortress_info.o fortress_util.o randlib.o fortress_random_t.o as63.o fortress_prior_t.o fortress_linalg.o filter.o fortress_model_t.o fortress_particles_t.o fortress_smc_particles_t.o fortress_particle_filter.o fortress_smc_t.o gensys.o

#LIBOBJS=fortress_prior_t.o fortress_model_t.o fortress_info.o randlib.o 

FC=mpif90 -O3 #-Wall -fcheck=all -g -fbacktrace
ifdef CONDA_BUILD
LIB=$(PREFIX)/lib
INC=$(PREFIX)/include
else
LIB=$(HOME)/anaconda3/lib
INC=$(HOME)/anaconda3/include
endif

#use export LD_LIBRARY_PATH=.
FPP=fypp
FRUIT=-I$(INC)/fruit -L$(LIB) -lfruit
FLAP=-I$(INC)/flap -L$(LIB) -lflap
#FORTRESS=-I/home/eherbst/Dropbox/code/fortress -L/home/eherbst/Dropbox/code/fortress -lfortress
FORTRESS=-I$(INC)/fortress -L$(LIB) -lfortress

JSON=-I$(INC)/json-fortran -L$(LIB)/json-fortran -ljsonfortran
.PHONY: all clean test test_library

%.o : %.f90
	$(FPP) -DGFORTRAN $< $(notdir $(basename $<))_tmp.f90
	$(FC) $(FRUIT) $(JSON) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) -o $(notdir $(basename $<)).o
	rm $(notdir $(basename $<))_tmp.f90

test_%.o : test_%.f90
	$(FPP) -DGFORTRAN $< $(notdir $(basename $<))_tmp.f90
	$(FC) $(FRUIT) $(JSON) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) $(FORTRESS) -o $(notdir $(basename $<)).o
	rm $(notdir $(basename $<))_tmp.f90


test_driver: test_driver.f90 $(LOBJS) test_model_t.o test_model.o test_prior.o test_random.o test_linalg.o test_smc.o test_util.o test_particles.o test_particle_filter.o test_gensys.o
	$(FC) $(FORTRESS) $^  -I. $(FRUIT) $(FLAP) $(FORTRESS) $(JSON) -lopenblas  -o $@ $(FORTRESS)

libfortress.so: fortress.f90  $(LIBOBJS) 
	$(FC) -shared -o $@  $^  

# test_model_t.o : test_model_t.f90 libfortress.so 
# 	$(FC) -c test/test_model_t.f90 -L. -lfortress -o test_model_t.o

# test_smc: smc_driver_mpi.f90 test_model_t. libfortress.so
# 	$(FC)  templates/smc_driver_mpi.f90 -L. test/test_model_t.f90 -lfortress $(FLAP) -llapack -o test_smc




test:
	python conda/run_test.py

clean:
	rm -f *.o *.mod
