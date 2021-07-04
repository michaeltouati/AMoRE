import rapc
import os

dir= os.path.dirname("figure/")
if not os.path.exists(dir):
    os.mkdir(dir)

# find the number of spatial grid points (Nx,Nz) by reading 'diag/nb[cm-3].dat'
t0  = []
x0  = []
nb_file0 = open('diag/nb[cm-3].dat', 'r')
line = nb_file0.readline()
line = line.strip()
array = line.split()
t0.append(float(array[0]))
x0.append(float(array[1]))
counter = 0
for line in nb_file0:
	line      = line.strip()
	array     = line.split()
	t0.append(float(array[0]))
	x0.append(float(array[1]))
	counter = counter + 1
	if x0[counter]!=x0[counter-1]:
		Nz = counter
		break
for line in nb_file0:
	line      = line.strip()
	array     = line.split()
	t0.append(float(array[0]))
	counter = counter + 1
	if t0[counter]!=t0[counter-1]:
		NxNz = counter
		break
nb_file0.close()	
Nx = NxNz / Nz

# find the number kinetic energy grid points Neps by reading 'diag/psi0_x[cm-3_keV-1].dat'
x0  = []
e0  = []
psi0_file0 = open('diag/psi0_x[cm-3_keV-1].dat', 'r')
line = psi0_file0.readline()
line = line.strip()
array = line.split()
x0.append(float(array[1]))
e0.append(float(array[2]))
counter = 0
for line in psi0_file0:
	line      = line.strip()
	array     = line.split()
	x0.append(float(array[1]))
	e0.append(float(array[2]))
	counter = counter + 1
	if x0[counter]!=x0[counter-1]:
		Neps = counter
		break
psi0_file0.close()	

dir= os.path.dirname("figure/distribution_function/")
if not os.path.exists(dir):
    os.mkdir(dir)

# Plot the distribution function at the location where the beam density is maximum 
# for a given transverse position xi
dir= os.path.dirname("figure/distribution_function/fb_x/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_distribution(Neps,e0,2,Nx,
								'diag/psi0_x[cm-3_keV-1].dat',
								'diag/psi1x_x[cm-3_keV-1].dat',
								'diag/psi1z_x[cm-3_keV-1].dat',
								'figure/distribution_function/fb_x/fb_x')
								
# Plot the distribution function at the location where the beam density is maximum 
# for a given depth zk
dir= os.path.dirname("figure/distribution_function/fb_z/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_distribution(Neps,e0,1,Nz,
								'diag/psi0_z[cm-3_keV-1].dat',
								'diag/psi1x_z[cm-3_keV-1].dat',
								'diag/psi1z_z[cm-3_keV-1].dat',
								'figure/distribution_function/fb_z/fb_z')
