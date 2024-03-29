#!/usr/bin/env bash
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

TEST_DIR=' Vlasov/HLL-1'
TEST_DIR+=' Vlasov/HLL-2'
TEST_DIR+=' Fokker-Planck/Implicit-collisions'
TEST_DIR+=' Parallelization/OpenMP'
TEST_DIR+=' MHD/Bi-temperature'
TEST_DIR+=' MHD/Magnetic-diffusion'
TEST_DIR+=' Laser-solid-interaction/Rear-side-refluxing'
TEST_DIR+=' Laser-solid-interaction/Both-sides-refluxing'
TEST_DIR+=' Solids/Al'
TEST_DIR+=' Solids/Cu'
TEST_DIR+=' Solids/Ta'
TEST_DIR+=' Solids/CH-like'
TEST_DIR+=' Solids/C-vitreous'
TEST_DIR+=' Plasmas/H'
TEST_DIR+=' Plasmas/Be'
TEST_DIR+=' Plasmas/Shock-ignition'
TEST_DIR+=' Spectrum/Gaussian'
TEST_DIR+=' Spectrum/Boltzmannian'
TEST_DIR+=' Spectrum/Bi-temperature'
TEST_DIR+=' Spectrum/Modified-bi-temperature'
TEST_DIR+=' Spectrum/Particle-In-Cell'
TEST_DIR+=' X-rays/Kalpha'
TEST_DIR+=' X-rays/Tracer'
# TEST_DIR+=' New-features/'

cd ../..
mv input-deck input-deck-old
for tst in ${TEST_DIR}; do \
	if [[ ${tst} = Solids/C-vitreous || ${tst} = Solids/CH-like ]]; then \
		mv sources/user/resistivity_tab.dat sources/user/resistivity_tab-old.dat ; \
		cp test-cases/Tests/${tst}/resistivity_tab.dat sources/user/ ; \
	fi
	if [ ${tst} = Spectrum/Particle-In-Cell ]; then \
		mv sources/user/spectrum_tab.dat sources/user/spectrum_tab-old.dat ; \
		cp test-cases/Tests/${tst}/spectrum_tab.dat sources/user/ ; \
	fi
	if [ ${tst} = Plasmas/Shock-ignition ]; then \
		mv sources/user/plasma_tab.dat sources/user/plasma_tab-old.dat ; \
		cp test-cases/Tests/${tst}/plasma_tab.dat sources/user/ ; \
	fi
	cp test-cases/Tests/${tst}/input-deck .
	./amore > output 
	mv output test-cases/Tests/${tst}/ ; \
	echo ${tst}'/output file generated'
	if [[ ${tst} = Solids/C-vitreous || ${tst} = Solids/CH-like ]]; then \
		mv sources/user/resistivity_tab-old.dat sources/user/resistivity_tab.dat ; \
	fi
	if [ ${tst} = Spectrum/Particle-In-Cell ]; then \
		mv sources/user/spectrum_tab-old.dat sources/user/spectrum_tab.dat ; \
	fi
	if [ ${tst} = Plasmas/Shock-ignition ]; then \
		mv sources/user/plasma_tab-old.dat sources/user/plasma_tab.dat ; \
	fi
done
rm -rf results/Vlasov output
mv input-deck-old input-deck
cd test-cases/Tests
