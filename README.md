## Written by Michaël J TOUATI

# AMoRE

Angular Momentum model Of Relativistic Electron beam transport (AMoRE) is a Fortran code parallelized using OpenMP for fast simulations of laser-generated relativistic electron beam transport in solids and dense plasmas. It computes the two first angular moments of the exact relativistic Belaiev-Budker kinetic equation for relativistic electron beam transport in solids and dense plasmas completed with the Minerbo maximum angular entropy closure. It thus takes into account collective effects with the self-generated electromagnetic fields as well as collisional effects with the slowing down of the electrons in collisions with plasmons, bound and free electrons and their angular scattering on both ions and electrons. This model allows for fast computations of relativistic electron beam transport while describing the kinetic distribution function evolution. More pieces of information can be found in my [PhD manuscript](https://tel.archives-ouvertes.fr/tel-01238782/document) and in this [peer-reviewed Physics paper](https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073014/pdf). Citations to these references are recommended and appreciated for publications of scientific results using AMoRE in peer-reviewed journals. This version of the code only implement 2D homogeneous cartesian spatial grids and OpenMP shared memory. 
The next release will implement 3D non-homogeneous cartesian meshes, MPI distributed memory and SIMD vectorization.

# Compiling the code

Modify the [makefile](https://github.com/michaeltouati/AMoRE/blob/main/Makefile) as a function of the wished compilation options and the compilers installed on your computer and then type :
```sh
make
```

# Running a simulation

Fill the wished input-deck [init_parameters](https://github.com/michaeltouati/AMoRE/blob/main/init_parameters) (all parameters are described inside) and type :
```sh
./amore
```
# Plotting the simulation results

Type for example :
```sh
python3 extract_tools/extract_maps.py
```
All simulation results are stored in text files located in the folder 'diags' and the python scripts used to make the plots are located in the folder 'extract_tools'.
The plots will be located in the folder 'figure'.
