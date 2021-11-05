[<img src='https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white' height="20">](https://fortran-lang.org/)
[<img src='https://img.shields.io/badge/GNU%20Make-black?style=for-the-badge&logo=gnu&logoColor=#7D929E' height="20">](https://www.gnu.org/software/make/)
[<img src='https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white' height="20">](https://www.gnu.org/software/bash/)
[<img src='https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54' height="20">](https://www.python.org/)
[<img src='https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white' height="20">](https://numpy.org/)
[<img src='https://matplotlib.org/_static/logo2_compressed.svg' height="20">](https://matplotlib.org/stable/index.html#)
[<img src='https://img.shields.io/badge/latex-%23008080.svg?style=for-the-badge&logo=latex&logoColor=white' height="20">](https://www.latex-project.org//)
[![GPLv3 License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![Downloads](https://img.shields.io/github/downloads/michaeltouati/AMoRE/total)](https://github.com/michaeltouati/AMoRE/releases)
[![Views : 14 days](https://img.shields.io/badge/dynamic/json?color=success&label=Views%20(<15%20days)&query=count&url=https://github.com/michaeltouati/AMoRE/blob/master/.github/view.json?raw=True&logo=github)](https://github.com/michaeltouati/AMoRE/actions/workflows/views.yml)
[![Clones : 14 days](https://img.shields.io/badge/dynamic/json?color=success&label=Clones%20(<15%20days)&query=count&url=https://github.com/michaeltouati/AMoRE/blob/master/.github/clone.json?raw=True&logo=github)](https://github.com/michaeltouati/AMoRE/actions/workflows/clones.yml)
[![Compilation check](https://github.com/michaeltouati/AMoRE/actions/workflows/compilation.yml/badge.svg?branch=master)](https://github.com/michaeltouati/AMoRE/actions/workflows/compilation.yml)
[![Tests check](https://github.com/michaeltouati/AMoRE/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/michaeltouati/AMoRE/actions/workflows/tests.yml)
[![Plotting tools check](https://github.com/michaeltouati/AMoRE/actions/workflows/plots.yml/badge.svg?branch=master)](https://github.com/michaeltouati/AMoRE/actions/workflows/plots.yml)
[![Open issues](https://img.shields.io/github/issues/michaeltouati/AMoRE)](https://github.com/michaeltouati/AMoRE/issues)
[![Closed issues](https://img.shields.io/github/issues-closed/michaeltouati/AMoRE)](https://github.com/michaeltouati/AMoRE/issues)
[![Open pull requests](https://img.shields.io/github/issues-pr/michaeltouati/AMoRE)](https://github.com/michaeltouati/AMoRE/pulls)
[![Closed pull requests](https://img.shields.io/github/issues-pr-closed/michaeltouati/AMoRE)](https://github.com/michaeltouati/AMoRE/pulls)
<!-- ![My Stats](https://github-readme-stats.vercel.app/api?username=michaeltouati&show_icons=true) -->

## Written by Michaël J TOUATI

<div style="text-align: justify">
[AMoRE](https://github.com/michaeltouati/AMoRE) (Angular Momentum model Of Relativistic Electron beam transport) is a Fortran code, parallelized with OpenMP and using Python3 for post-processing, that allows for the fast simulation of a relativistic electron beam transport in a denser solid or plasma while describing the full beam electron phase-space evolution. 

To achieve this, it computes the two first momentum angular moments of the beam electrons Vlasov-Fokker-Planck-Belaiev-Budker (VFPBB) kinetic equation (the two first order equations of the kinetic equation Cartesian tensor scalar product expansion), completed with the Minerbo maximum angular entropy closure to express the distribution function second order momentum angular moment needed in the first order equation. The resulting kinetic equations are coupled with the target electrons and ions Magneto-Hydrodynamic (MHD) equations considering time scales on which the target electron return current has already set up. [AMoRE](https://github.com/michaeltouati/AMoRE) thus takes into account both collective effects, with the self-generated electromagnetic fields, and collisional effects, with the slowing down of beam electrons in collisions with plasmons, bound and free electrons and their angular scattering on both ions and electrons. The kinetic energy cut-off separating target electrons and beam electrons is assumed to be around 10 keV for the beam electrons collisional friction and diffusion VFPBB terms and the solid/plasma species MHD equations to be valid. More pieces of information can be found in my [PhD manuscript](https://tel.archives-ouvertes.fr/tel-01238782/document) and in this [peer-reviewed article](https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073014/pdf). In the latter, there is no mention of updates in relation with the possible beam electrons specular reflection boundary conditions (refluxing), thermal capacities valid at ambient temperatures according to the Einstein and Debye models as well as the distinction between collisions of free electrons with s, p or d-band electrons in the electron-electron collisions contribution on transport coefficients. Citations to these references are recommended and appreciated for publications of scientific results using [AMoRE](https://github.com/michaeltouati/AMoRE) in peer-reviewed journals. 

Python scripts, using the Matplotlib and Numpy packages, are provided to automatically extract and plot the simulation results. The simulation parameters are described in the [input-deck](https://github.com/michaeltouati/AMoRE/blob/main/input-deck) and they can be modified without having to recompile the code. Compilation rules can be modified in the [makefile](https://github.com/michaeltouati/AMoRE/blob/main/makefile) depending on the user compiler preferences. Tools for testing the compilation of the code and tools for checking the simulation parameters are provided. X-ray Kα and Kβ photon emission from the target can be computed using a provided [X-ray Kα and Kβ table](https://github.com/michaeltouati/AMoRE/blob/master/sources/data/Kalpha_tab.dat). The user can use its own [tabulated electrical resistivity](https://github.com/michaeltouati/AMoRE/blob/master/sources/user/resistivity_tab.dat) for the solid/plasma and/or its own [tabulated plasma material-density-temperature profile](https://github.com/michaeltouati/AMoRE/blob/master/sources/user/plasma_tab.dat) and/or its own [tabulated electron beam kinetic energy spectrum](https://github.com/michaeltouati/AMoRE/blob/master/sources/user/spectrum_tab.dat).
<div>
  
# Compiling the code

Modify the [makefile](https://github.com/michaeltouati/AMoRE/blob/main/makefile) as a function of the Fortran compiler installed on your computer and then type
```sh
make
```
The compilation can be tested by typing
```sh
make test
```
The tests consist in comparing file1 and file2 where :
* file1 is one test simulation terminal output performed with an input deck located in the directory 'test-cases/Tests/' and
* file2 is the terminal output of the corresponding simulation already performed by the developper also located in 'test-cases/Tests/'.

# Running a simulation

Fill the wished [input-deck](https://github.com/michaeltouati/AMoRE/blob/main/input-deck) (all parameters are described inside), eventually check them by typing
```sh
./check-input-deck
```
or
```sh
make check
```
and then type
```sh
./amore
```
or
```sh
make run
```

# Plotting the simulation results

Type
```sh
make plot
```
to plot all the simulation results that are stored in the directory 'results'. It can also be plotted separately all 2D spatial density plots of relativistic electron beam and target species hydrodynamic moments as well as the self-generated electromagnetic fields and target transport coefficients at each damped time step by typing
```sh
make plot_maps
```
or
```sh
python3 sources/plot/plot_maps.py
```
Similarly, the relativistic beam electron 2D-3P phase-space distribution can be plotted at the maximum beam density location and at each damped time step by typing
```sh
make plot_distribution
```
The simulation can be checked by typing
```sh
make plot_initialization
```
```sh
make plot_material
```
and
```sh
make plot_energy
```
that plots respectively :
* the initialiazed properties of the relativistic electron beam (longitudinal and transversal spatial distributions, momentum angular distribution as a function of transverse locations and kinetic energy spectrum),
* the transport coefficients of the solid/liquid/Warm-Dense-Matter/plasma where the beam is propagating through (ionization state, thermal capacities, electrical and thermal conductivities as well as electron-ion/lattice coupling factor) and 
* the time evolution of the different energy contributions in the beam transport (injected electron beam energy, collisional and collective electron beam energy losses, instantaneous electron beam energy in the target and electron beam energy escaping from the target). 

The python scripts are located in the directory [sources/plot](https://github.com/michaeltouati/AMoRE/tree/main/sources/plot).
The plots will be located in the directory 'figures'.

# License
[AMoRE](https://github.com/michaeltouati/AMoRE) is distributed under the terms of the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license. 
