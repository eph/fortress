SRC=src
TEST=test
VPATH=.:$(SRC):$(TEST):templates

LIBOBJS=fortress_info.o fortress_util.o randlib.o fortress_random_t.o as63.o fortress_prior_t.o fortress_linalg.o filter.o  fortress_model_t.o fortress_particles_t.o fortress_smc_t.o 

#LIBOBJS=fortress_prior_t.o fortress_model_t.o fortress_info.o randlib.o 

FC=mpif90 -fbounds-check
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

.PHONY: all clean test test_library

%.o : %.f90
	$(FPP) -DGFORTRAN $< $(notdir $(basename $<))_tmp.f90
	$(FC) $(FRUIT) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) -o $(notdir $(basename $<)).o
	rm $(notdir $(basename $<))_tmp.f90

test_prior.o : test_prior.f90
	$(FC) $(FRUIT) -c $<

test_random.o : test_random.f90
	$(FC) $(FRUIT) -c $<

test_library: test_library.f90 libfortress.so
	$(FC) src/test_library.f90  -L. -lfortress -llapack $(FLAP) -o test_library

test_driver: test_driver.f90 $(LOBJS) test_model_t.o test_model.o test_prior.o test_random.o test_linalg.o test_smc.o test_util.o test_particles.o
	$(FC) $^  -I. $(FRUIT) $(FLAP) -L. -lfortress  -llapack  -o $@ 

libfortress.so: fortress.f90  $(LIBOBJS) 
	$(FC) -shared -o $@  $^  

#test_model_t.o : test_model_t.f90 libfortress.so 
# 	$(FC) -c test/test_model_t.f90 -L. -lfortress -o test_model_t.o

test_smc: smc_driver_mpi.f90 test_model_t.f90 libfortress.so
	$(FC)  templates/smc_driver_mpi.f90 -L. test/test_model_t.f90 -lfortress $(FLAP) -llapack -o test_smc 



test:
	python conda/run_test.py

clean:
	rm -f *.o *.mod
