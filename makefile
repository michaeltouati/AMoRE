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

######################
######################
##  INTEL compiler  ##
######################
######################

# F90 = ifort

##########
# openMP #
##########

# OPTS = -O3 -fopenmp -r8
# OPTS = -O2 -fopenmp

#########
# debug #
#########

# OPTS = -O -g -debug all -debug-parameters all -C -fopenmp
# OPTS = -g -traceback -fopenmp -r8 -std95 -fpe0 -debug all -debug-parameters all -C

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
# OPTS = -fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow

#####################################
#####################################
##        Compilation rules        ##
#####################################
#####################################

SRC_PATH    = sources/
SRC_PATH_PY = sources/plot/

SRCS = acuracy.f90 \
	constants.f90 \
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

#############
#############
## Running ##
#############
#############

run :
	@./amore

##############
##############
## Cleaning ##
##############
##############

clean_results :
	@rm -rf results

clean_figures :
	@rm -rf figures

clean_dist :
	@rm -rf *.mod *.o amore ${SRC_PATH_PY}__pycache__

clean : clean_results clean_figures clean_dist

##############
##############
## Plotting ##
##############
##############

plot_energy :
	@python3 ${SRC_PATH_PY}plot_energy.py

plot_initialization :
	@python3 ${SRC_PATH_PY}plot_initialization.py

plot_material :
	@python3 ${SRC_PATH_PY}plot_material.py

plot_maps :
	@python3 ${SRC_PATH_PY}plot_maps.py

plot_distribution :
	@python3 ${SRC_PATH_PY}plot_distribution.py

plot : plot_energy \
	plot_initialization \
	plot_material \
	plot_maps \
	plot_distribution

###############
###############
##  Testing  ##
###############
###############

TEST_DIR  = Vlasov/HLL-1
TEST_DIR += Vlasov/HLL-2
TEST_DIR += Parallelization/OpenMP
TEST_DIR += Fokker-Planck/Implicit-collisions
# TEST_DIR += Boundary-conditions/Periodic
# TEST_DIR += Boundary-conditions/Absorbing
# TEST_DIR += Vlasov-linear/Donor-cell
# TEST_DIR += Vlasov-linear/Lax-Wendroff
# TEST_DIR += Vlasov-linear/Beam-Warming
# TEST_DIR += Vlasov-linear/Fromm
# TEST_DIR += Vlasov-nonlinear/Minmod
# TEST_DIR += Vlasov-nonlinear/Superbee
# TEST_DIR += Vlasov-nonlinear/Van-Leer
# TEST_DIR += Vlasov-nonlinear/MUSCL1
# TEST_DIR += Vlasov-nonlinear/MUSCL2
# TEST_DIR += Academic-cases/Landau
# TEST_DIR += Academic-cases/Wakefield
# TEST_DIR += Academic-cases/Two-stream-insta
# TEST_DIR += New-features/Example

RED   =$(shell tput setaf 1)
GREEN =$(shell tput setaf 2)
RESET =$(shell tput sgr0)
TESTS:=$(sort ${TEST_DIR})
test :
	@mv input-deck input-deck-old
	@echo '---------------------------------'
	@echo '        TESTS DESCRIPTION        '
	@echo '---------------------------------'
	@echo '                                 '
	@echo 'The tests consist in performing  '
	@echo ' diff file1 file2 where :        '
	@echo ' * file1 is the test simulation  '
	@echo '   terminal output               '
	@echo ' * file2 the terminal output of  '
	@echo '   the corresponding simulation  '
	@echo '   performed by the developper   '
	@echo '   located in test-cases/Tests/  '
	@echo '                                 '
	@for tst in ${TESTS}; do \
		echo '---------------------------------' ; \
		echo $${tst}' :'  ; \
		cp test-cases/Tests/$${tst}/input-deck . ; \
		./amore > test.output ; \
		if hash tac 2>/dev/null; then \
			tail -n +0 test.output | head -n -4 > file1 ; \
			tail -n +0 test-cases/Tests/$${tst}/output | head -n -4 > file2 ; \
		else \
			tail -n +0 test.output | tail -r | tail -n +5 | tail -r > file1 ; \
			tail -n +0 test-cases/Tests/$${tst}/output | tail -r | tail -n +5 | tail -r > file2 ; \
		fi ; \
		diff file1 file2; \
		TST=$$?; \
		rm file1; rm file2; \
		if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; \
	done
	@rm -rf results/Vlasov
	@mv input-deck-old input-deck
	@rm -f  test.output