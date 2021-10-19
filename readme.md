## Written by Michaël J TOUATI

[AMoRE](https://github.com/michaeltouati/AMoRE) (Angular Momentum model Of Relativistic Electron beam transport) is a unique Fortran code parallelized using OpenMP allowing for fast simulations of relativistic electron beam transport in solid or dense plasma targets while describing the full relativistic beam electron phase-space evolution. It computes the two first angular moments of the relativistic beam electrons Vlasov-Fokker-Planck-Belaiev-Budker kinetic equation, completed with the Minerbo maximum angular entropy closure, that is coupled with the target electrons and ions Magneto-Hydrodynamic equations considering time scales on which the return current has already set up. It thus takes into account both collective effects, with the self-generated electromagnetic fields, and collisional effects, with the slowing down of beam electrons in collisions with plasmons, bound and free electrons and their angular scattering on both ions and electrons. More pieces of information can be found in my [PhD manuscript](https://tel.archives-ouvertes.fr/tel-01238782/document) and in this [peer-reviewed article](https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073014/pdf). Citations to these references are appreciated for publications of scientific results using AMoRE in peer-reviewed journals. 

# Compiling the code

Modify the [makefile](https://github.com/michaeltouati/AMoRE/blob/main/Makefile) as a function of the Fortran compiler installed on your computer and then type
```sh
make
```

# Running a simulation

Fill the input-deck [init_parameters](https://github.com/michaeltouati/AMoRE/blob/main/init_parameters) (all parameters are described inside) and type :
```sh
./amore
```
# Plotting the simulation results

Type
```sh
make plot
```
to plot all the simulation results that are stored in the directory 'diag'. It can also be plotted separately all 2D spatial density plots of relativistic electron beam and target species hydrodynamic moments as well as the self-generated electromagnetic fields and target transport coefficients at each damped time step by typing
```sh
python3 sources/tools/plot_maps.py
```
Similarly, the relativistic beam electron 2D-3P phase-space distribution can be plotted at the maximum beam density location and at each damped time step by typing
```sh
python3 sources/tools/plot_distribution.py
```
The simulation can be checked by typing 
```sh
python3 sources/tools/plot_initialization.py
```
```sh
python3 sources/tools/plot_material.py
```
and
```sh
python3 sources/tools/plot_energy.py
```
that plots respectively :
* the initialiazed properties of the relativistic electron beam (longitudinal and transversal spatial distributions, momentum angular distribution as a function of transverse locations and kinetic energy spectrum),
* the transport coefficients of the solid/liquid/Warm-Dense-Matter/plasma where the beam is propagating through (ionization state, thermal capacities, electrical and thermal conductivities as well as electron-ion/lattice coupling factor) and 
* the time evolution of the different energy contributions in the beam transport (injected electron beam energy, collisional and collective electron beam energy losses, instantaneous electron beam energy in the target and electron beam energy escaping from the target). 

The python scripts are located in the directory [sources/tools](https://github.com/michaeltouati/AMoRE/tree/main/sources/tools).
The plots will be located in the directory 'figures'.
