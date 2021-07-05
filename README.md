# AMoRE
Dr M. Touati - UCLA - 2015/11/28 - currently at the CLPU - mtouati@clpu.es

Angular Momentum model Of Relativistic Electron beam transport in solids or dense plasmas (AMoRE) is a hybrid reduced model for laser-generated relativistic electron beam transport in solids and dense plasmas. It is based on the two first angular moments of the relativistic Belyaev-Budker kinetic equation completed with the Minerbo maximum angular entropy closure. It takes into account collective effects with the self-generated electromagnetic fields as well as collisional effects with the slowing down of the electrons in collisions with plasmons, bound and free electrons and their angular scattering on both ions and electrons. This model allows for fast computations of relativistic electron beam transport while describing the kinetic distribution function evolution. More pieces of information can be found here :
https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073014/pdf
and 
https://tel.archives-ouvertes.fr/tel-01238782/document

For publications of scientific results using AMoRE in peer-review journals, do not forget to cite the latter article and/or PhD manuscript. 

This version of the code only implement 2D homogeneous cartesian spatial grids and OpenMP shared memory. 
The next release will implement 3D non-homogeneous cartesian meshes, MPI distributed memory and SIMD vectorization.

In order to :

1) Compile the code :

Modify the makefile as a function of the wished compilation options and the compilers installed on your computer and then type :

make

2) Run a simulation :

Fill the wished input-deck 'init_parameters' (all parameters are described inside) and type :

./amore

3) Plot the simulation results :

Type for example :

python3 extract_tools/extract_maps.py

All simulation results are stored in text files located in the folder 'diags' and the python scripts used to make the plots are located in the folder 'extract_tools'.
The plots will be located in the folder 'figure'.
