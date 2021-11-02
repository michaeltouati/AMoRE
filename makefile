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
# OPTS = -fopenmp 

#########
# debug #
#########

# OPTS = -O -g -debug all -debug-parameters all -C -fopenmp
# OPTS = -O0 -traceback -r8 -std08 -fpe0 -g -debug all -debug-parameters all -C -fopenmp

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

# OPTS = -O3 -fopenmp -fdefault-real-8
OPTS = -O2 -fopenmp

#########
# debug #
#########

#OPTS = -O3 -Wall -fcheck=all -g -fbacktrace -fopenmp
# OPTS = -fdefault-real-8 -O0 -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow

#####################################
#####################################
##        Compilation rules        ##
#####################################
#####################################

SRC_PATH    = sources/
SRC_PATH_PY = sources/plot/

SRCS = Acuracy.f90 \
	Constants.f90 \
	Physics_library.f90 \
	Fokker_Planck_coef.f90 \
	Input.f90 \
	Transport_coef.f90 \
	VFP.f90 \
	Initialization.f90 \
	Diagnostics.f90 \
	EM_fields.f90 \
	Heat_equations.f90 \
	AMoRE.f90

SRCS_CHK = Acuracy.f90 \
	Constants.f90 \
	Physics_library.f90 \
	Fokker_Planck_coef.f90 \
	Input.f90 \
	Input_chk.f90

OBJTS := $(SRCS:%.f90=%.o)

OBJTS_CHK := $(SRCS_CHK:%.f90=%.o)

all : check-input-deck amore \
	remove-compilation-files

amore : $(OBJTS)
	$(F90) $(OPTS) $(OBJTS) -o amore

%.o : $(SRC_PATH)%.f90
	$(F90) $(OPTS) -c $(SRC_PATH)$*.f90

check-input-deck : $(OBJTS_CHK)
	$(F90) $(OPTS) $(OBJTS_CHK) -o check-input-deck

remove-compilation-files :
	@rm *.mod
	@rm *.o

#########################
#########################
## Checking input-deck ##
#########################
#########################

check :
	@./check-input-deck

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
	@rm -rf *.mod *.o amore check-input-deck ${SRC_PATH_PY}__pycache__

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
TEST_DIR += MHD/Bi-temperature
TEST_DIR += MHD/Magnetic-diffusion
TEST_DIR += Laser-solid-interaction/Rear-side-refluxing
TEST_DIR += Laser-solid-interaction/Both-sides-refluxing
TEST_DIR += Solids/Al
TEST_DIR += Solids/Cu
TEST_DIR += Solids/Ta
TEST_DIR += Solids/CH-like
TEST_DIR += Solids/C-vitreous
TEST_DIR += Plasmas/H
TEST_DIR += Plasmas/Be
# TEST_DIR += Plasmas/Tabulated
TEST_DIR += Spectrum/Gaussian
TEST_DIR += Spectrum/Boltzmannian
TEST_DIR += Spectrum/Bi-temperature
TEST_DIR += Spectrum/Modified-bi-temperature
# TEST_DIR += Spectrum/Tabulated
TEST_DIR += X-rays/Kalpha
TEST_DIR += X-rays/Tracer
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
		if [ $${tst} = Solids/C-vitreous ]; then \
			mv sources/user/resistivity_tab.dat sources/user/resistivity_tab-old.dat ; \
			cp test-cases/Tests/Solids/C-vitreous/resistivity_tab.dat sources/user/ ; \
		fi ; \
		if [ $${tst} = Solids/CH-like ]; then \
			mv sources/user/resistivity_tab.dat sources/user/resistivity_tab-old.dat ; \
			cp test-cases/Tests/Solids/CH-like/resistivity_tab.dat sources/user/ ; \
		fi ; \
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
		if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; exit 2; fi; \
		if [ $${tst} = Solids/C-vitreous ]; then \
			mv sources/user/resistivity_tab-old.dat sources/user/resistivity_tab.dat ; \
		fi ; \
		if [ $${tst} = Solids/CH-like ]; then \
			mv sources/user/resistivity_tab-old.dat sources/user/resistivity_tab.dat ; \
		fi ; \
	done
	@rm -rf results/Vlasov
	@mv input-deck-old input-deck
	@rm -f  test.output