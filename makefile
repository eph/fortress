SRC=src
TEST=test
VPATH=.:$(SRC):$(TEST):$(SRC)/external

LIBOBJS=fortress_info.o fortress_util.o fortress_linalg.o randlib.o fortress_random_t.o as63.o fortress_prior_t.o filter.o fortress_model_t.o fortress_VAR_t.o fortress_particles_t.o fortress_smc_particles_t.o fortress_particle_filter.o fortress_smc_t.o gensys.o


FC = gfortran
ifeq ($(FC), ifort)
	FC=mpif90 -mkl -openmp -O3 -nocheck -inline-level=2 -shared-intel -mcmodel=medium -xSSE4.2 -ipo
	FCDEC=-DIFORT 
endif

ifeq ($(FC), gfortran)
#	FC=mpif90 -f90=gfortran -ffree-line-length-1000 -Wall -fcheck=all -g -fbacktrace #03
	FC=mpif90 -f90=gfortran -O3 -ffast-math -ffree-line-length-1000 #-Wall -fcheck=all -g -fbacktrace #03
	FCDEC=-DGFORTRAN
endif



ifdef CONDA_BUILD
LIB=$(PREFIX)/lib
INC=$(PREFIX)/include
else
LIB=$(HOME)/anaconda3/lib
INC=$(HOME)/anaconda3/include
endif

#use export LD_LIBRARY_PATH=.
FPP=fypp
FRUIT=-I$(INC)/fruit -L$(LIB) -lfruit -Wl,-rpath=$(LIB)
FLAP=-I$(INC)/flap -L$(LIB) -lflap
#FORTRESS=-I$(INC)/fortress -L$(LIB) -lfortress
FORTRESS= -L. -lfortress -Wl,-rpath=.
JSON=-I$(INC)/json-fortran -L$(LIB)/json-fortran -ljsonfortran


.PHONY: all clean test test_library

%.o : %.f90
	$(FPP) $(FCDEC) $< $(notdir $(basename $<))_tmp.f90
	$(FC) $(FRUIT) $(JSON) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) -o $(notdir $(basename $<)).o
	rm $(notdir $(basename $<))_tmp.f90

test_%.o : test_%.f90
	$(FPP) $(FCDEC) $< $(notdir $(basename $<))_tmp.f90
	$(FC)  $(FORTRESS) $(FRUIT) $(JSON) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) -o $(notdir $(basename $<)).o 
	rm $(notdir $(basename $<))_tmp.f90


test_driver: test_driver.f90 $(LOBJS) test_model_t.o test_VAR_t.o test_VAR.o test_model.o test_prior.o test_random.o test_linalg.o test_smc.o test_util.o test_particles.o test_particle_filter.o test_gensys.o
	$(FC) $(FORTRESS) $(FLAP) $(FRUIT) $(JSON) $^ -o $@ -llapack 

libfortress.so: fortress.f90  $(LIBOBJS) 
	$(FC) -shared -o $@  $^ 

test: test_driver.o
	python conda/run_test.py

clean:
	rm -f *.o *.mod
