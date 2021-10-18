#######################################################################
##                                                                   ##
## Angular Momentum Model Of Relativistic Electron beam (AMoRE) code ##
##                                                                   ##
## Copyright © 2015 Michaël J TOUATI                                 ##
##                                                                   ##
## This file is part of AMoRE.                                       ##
##                                                                   ##
## AMoRE is free software: you can redistribute it and/or modify     ##
## it under the terms of the GNU General Public License as published ##
## by the Free Software Foundation, either version 3 of the License, ##
## or (at your option) any later version.                            ##
##                                                                   ##
## AMoRE is distributed in the hope that it will be useful,          ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of    ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
## GNU General Public License for more details.                      ##
##                                                                   ##
## You should have received a copy of the GNU General Public License ##
## along with AMoRE. If not, see <https://www.gnu.org/licenses/>.    ##
##                                                                   ##
#######################################################################
## Initial commit written by Michaël J TOUATI - Oct. 2015

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

# F90 = mpif90
F90 = gfortran

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

SRC_PATH    = sources/
SRC_PATH_PY = sources/tools/

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
	AMoRE.f90

OBJTS := $(SRCS:%.f90=%.o)

all : amore

amore : $(OBJTS)
	$(F90) $(OPTS) $(OBJTS) -o amore
	rm *.mod
	rm *.o

%.o : $(SRC_PATH)%.f90
	$(F90) $(OPTS) -c $(SRC_PATH)$*.f90

clean_results :
	@rm -rf diag/*

clean_figures :
	@rm -rf figures/*

clean : clean_results clean_figures
	@rm -rf *.mod *.o sources/*.o sources/*.mod amore ${SRC_PATH_PY}__pycache__

plot :
	@python3 ${SRC_PATH_PY}extract_energy.py
	@python3 ${SRC_PATH_PY}extract_initialization.py
	@python3 ${SRC_PATH_PY}extract_material.py
	@python3 ${SRC_PATH_PY}extract_maps.py
	@python3 ${SRC_PATH_PY}extract_distribution.py