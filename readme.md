## Written by Michaël J TOUATI

# AMoRE

Angular Momentum model Of Relativistic Electron beam transport (AMoRE) is a Fortran code parallelized using OpenMP for fast simulations of laser-generated relativistic electron beam transport in solids or dense plasmas. It computes the two first angular moments of the exact relativistic Vlasov-Fokker-Planck-Belaiev-Budker kinetic equation completed with the Minerbo maximum angular entropy closure and coupled to Magneto-Hydrodynamic equations considering time scales during which the return current has already set up. It thus takes into account collective effects with the self-generated electromagnetic fields as well as collisional effects with the slowing down of beam electrons in collisions with plasmons, bound and free electrons and their angular scattering on both ions and electrons. This model allows for fast computations while also describing the relativistic beam electron phase-space evolution. More pieces of information can be found in my [PhD manuscript](https://tel.archives-ouvertes.fr/tel-01238782/document) and in this [peer-reviewed article](https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073014/pdf). Citations to these references are recommended and appreciated for publications of scientific results using AMoRE in peer-reviewed journals. 

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
to plot all the simulation results automatically. It can also be plotted separately
all 2D density plots of both relativistic beam electrons and target electrons/ions hydrodynamic moments as well as the electromagnetic fields and target transport coefficient by typing
```sh
python3 sources/tools/extract_maps.py
```
Similarly, the relativistic beam electron 2D-3P phase-space distribution can be plotted at the maximum beam density location by typing
```sh
python3 sources/tools/extract_distribution.py
```
and the simulation can be checked by typing 
```sh
python3 sources/tools/extract_initialization.py
```
```sh
python3 sources/tools/extract_material.py
```
and
```sh
python3 sources/tools/extract_energy.py
```
that plot respectively :
* the initial properties of the relativistic electron beam (longitudinal and transversal spatial distribution, momentum angular distribution and spectrum),
* the transport coefficients of the medium where the beam is propagating through (ionization state, thermal capacities, electrical and thermal conductivities as well as electron-ion/lattice coupling factor) and 
* the time evolution of the different energy contribution in the simulation (injected electron beam energy, collisional and collective electron beam energy losses, instantaneous electron beam energy in the medium and escaping electron beam energy). 

The python scripts used to make the plots are located in the directory [sources/tools](https://github.com/michaeltouati/AMoRE/tree/main/sources/tools).
The plots will be located in the folder 'figure'.
