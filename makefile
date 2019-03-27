SRC=src
TEST=test
VPATH=.:$(SRC):$(TEST):$(SRC)/external

LIBOBJS=fortress_info.o fortress_constants.o fortress_util.o fortress_linalg.o randlib.o fortress_random_t.o as63.o fortress_prior_t.o filter.o fortress_model_t.o fortress_VAR_t.o fortress_particles_t.o fortress_smc_particles_t.o fortress_particle_filter.o fortress_smc_utilities.o fortress_smc_t.o gensys.o 


# ifeq ($(FC), ifort)
# 	FC=mpif90 -mkl -openmp -O3 -nocheck -inline-level=2 -shared-intel -mcmodel=medium -xSSE4.2 -ipo
# 	FCDEC=-DIFORT 
# endif


COMPILER=mpif90 -fstack-protector -fimplicit-none -O3 -ffree-line-length-1000 # -Wall -fcheck=all -g -fbacktrace #03
#COMPILER=mpif90 -Wall -fcheck=all -g -fbacktrace -ffree-line-length-1000 # 
FCDEC=-DGFORTRAN




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
	$(COMPILER) $(FRUIT) $(JSON) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) -o $(notdir $(basename $<)).o
	rm $(notdir $(basename $<))_tmp.f90

test_%.o : test_%.f90
	$(FPP) $(FCDEC) $< $(notdir $(basename $<))_tmp.f90
	$(COMPILER)  $(FORTRESS) $(FRUIT) $(JSON) -fPIC -c $(notdir $(basename $<)_tmp.f90) $(FLAP) -o $(notdir $(basename $<)).o 
	rm $(notdir $(basename $<))_tmp.f90


test_driver: test_driver.f90 libfortress.so $(LOBJS) test_model_t.o test_model.o test_prior.o test_random.o test_linalg.o test_smc.o test_util.o test_particles.o test_particle_filter.o test_gensys.o test_json.o test_model_circle_t.o
	$(COMPILER) $^ -o $@ $(JSON) $(FORTRESS) -lopenblas $(FLAP) $(JSON) $(FRUIT) 

#smc_driver: smc_driver.f90 /home/eherbst/Dropbox/var_smc_estimation/replication-code/smc_msvar/_fortress_tmp/model_t.f90 $(LIBOBJS)
smc_driver: smc_driver.f90 test/test_model_circle_t.f90 $(LIBOBJS)
	$(COMPILER) $(FORTRESS) $(FLAP) $(FRUIT) $(JSON) $^ -o $@ -llapack $(FORTRESS) $(FLAP) $(FRUIT) $(JSON) -lopenblas $(FORTRESS)


libfortress.so: fortress.f90  $(LIBOBJS) 
	$(COMPILER) -shared -o $@  $^ 


test_function_enclosure : test_function_enclosure.f90 libfortress.so test_model_t.o
	$(COMPILER) $^ -o $@ $(JSON) $(FORTRESS) -lopenblas $(FLAP) $(JSON) $(FRUIT)
test: test_driver.o
	python conda/run_test.py

clean:
	rm -f *.o *.mod
