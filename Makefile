#############
#############
##  ifort  ##
#############
#############

#F90 = ifort

##########
# openMP #
##########

#OPTS = -O3 -openmp -openmp-report2 -r8

#########
# debug #
#########

#OPTS = -O -g -debug all -debug-parameters all\
#	 -C -openmp -openmp-report2

################
################
##  gfortran  ##
################
################

F90 = mpif90

##########
# openMP #
##########

#OPTS = -O3 -fopenmp -fdefault-real-8
OPTS = -O2 -fopenmp

#########
# debug #
#########

#OPTS = -O3 -Wall -fcheck=all -g -fbacktrace -fopenmp

#####################################
#####################################
##        Compilation rules        ##
#####################################
#####################################

SRC_PATH    = source/
SRC_PATH_PY = source/extract/

SRCS = constants.f90 \
	physics_library.f90 \
	Fokker_Planck_coef.f90 \
	input.f90 \
	transport_coef.f90 \
	VFP.f90 \
	initialization.f90 \
	diagnostics.f90 \
	EM_fields.f90 \
	temperatures.f90 \
	main.f90

OBJTS := $(SRCS:%.f90=%.o)

all : AMoRE

AMoRE : $(OBJTS)
	$(F90) $(OPTS) $(OBJTS) -o AMoRE
	rm *.mod
	rm *.o

%.o : $(SRC_PATH)%.f90
	$(F90) $(OPTS) -c $(SRC_PATH)$*.f90

clean :
	rm -rf *.mod *.o source/*.o source/*.mod diag/*.dat figure/* AMoRE

clean_results :
	rm -rf diag/*.dat

clean_figures :
	rm -rf figure/*
