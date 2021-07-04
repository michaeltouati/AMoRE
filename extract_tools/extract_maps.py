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

# plot 2D maps
dir= os.path.dirname("figure/2D_maps/")
if not os.path.exists(dir):
    os.mkdir(dir)

print('Target ion density ni')
dir= os.path.dirname("figure/2D_maps/ni/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/ni[cm-3].dat',
						Nx,Nz,
						'jet',
						r'$\mathbf{n_i}\,(\mathrm{cm}^{-3})$',
						'figure/2D_maps/ni/ni_',
						0)

print('Target electron density ne')
dir= os.path.dirname("figure/2D_maps/ne/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/ne[cm-3].dat',
						Nx,Nz,
						'jet',
						r'$\mathbf{n_e}\,(\mathrm{cm}^{-3})$',
						'figure/2D_maps/ne/ne_',
						0)

print('Beam density nb')
dir= os.path.dirname("figure/2D_maps/nb/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/nb[cm-3].dat',
						Nx,Nz,
						'Blues',
						r'$\mathbf{n_b}\,(\mathrm{cm}^{-3})$',
						'figure/2D_maps/nb/nb_',
						0)

print('Beam current density jbz')
dir= os.path.dirname("figure/2D_maps/jbz/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/jb_z[A_cm-2].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{j_{b,z}}\,(\mathrm{A/cm}^{2})$',
						'figure/2D_maps/jbz/jbz_',
						0)

print('Beam current density jbx') 
dir= os.path.dirname("figure/2D_maps/jbx/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/jb_x[A_cm-2].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{j_{b,x}}\,(\mathrm{A/cm}^{2})$',
						'figure/2D_maps/jbx/jbx_',
						0)
 
print('Beam current density |jb|')    
dir= os.path.dirname("figure/2D_maps/jb/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh_abs('diag/jb_z[A_cm-2].dat',
								'diag/jb_x[A_cm-2].dat',
								Nx,Nz,
								'Reds',
								r'$|\mathbf{j_b}\,(\mathrm{A/cm}^{2})|$',
								'figure/2D_maps/jb/jb_')

print('Return current density jez')
dir= os.path.dirname("figure/2D_maps/jez/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/je_z[A_cm-2].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{j_{e,z}}\,(\mathrm{A/cm}^{2})$',
						'figure/2D_maps/jez/jez_',
						0)

print('Return current density jex') 
dir= os.path.dirname("figure/2D_maps/jex/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/je_x[A_cm-2].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{j_{e,x}}\,(\mathrm{A/cm}^{2})$',
						'figure/2D_maps/jex/jex_',
						0)
 
print('Return current density |je|')    
dir= os.path.dirname("figure/2D_maps/je/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh_abs('diag/je_z[A_cm-2].dat',
								'diag/je_x[A_cm-2].dat',
								Nx,Nz,
								'Reds',
								r'$|\mathbf{j_e}\,(\mathrm{A/cm}^{2})|$',
								'figure/2D_maps/je/je_')

print('Electric field Ex')
dir= os.path.dirname("figure/2D_maps/Ex/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/E_x[V_m-1].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{E_{x}}\,(\mathrm{V/m})$',
						'figure/2D_maps/Ex/Ex_',
						0)

print('Electric field Ez')
dir= os.path.dirname("figure/2D_maps/Ez/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/E_z[V_m-1].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{E_{z}}\,(\mathrm{V/m})$',
						'figure/2D_maps/Ez/Ez_',
						0)

print('Electric field |E|')    
dir= os.path.dirname("figure/2D_maps/E/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh_abs('diag/E_z[V_m-1].dat',
								'diag/E_x[V_m-1].dat',
								Nx,Nz,
								'Reds',
								r'$|\mathbf{E}\,(\mathrm{V/m})|$',
								'figure/2D_maps/E/E_')

print('Magnetic field By')
dir= os.path.dirname("figure/2D_maps/By/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/B_y[Tesla].dat',
						Nx,Nz,
						'seismic',
						r'$\mathbf{B_y}\,(\mathrm{T})$',
						'figure/2D_maps/By/By_',
						0)

print('Density of power deposited on background electrons We')
dir= os.path.dirname("figure/2D_maps/We/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/We[erg_s-1_cm-3].dat',
						Nx,Nz,
						'OrRd',
						r'$\mathbf{W_e}\,(\mathrm{erg/s.cm}^{3})$',
						'figure/2D_maps/We/We_',
						0)
						
print('Density of power deposited on background ions Wi')
dir= os.path.dirname("figure/2D_maps/Wi/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/Wi[erg_s-1_cm-3].dat',
						Nx,Nz,
						'OrRd',
						r'$\mathbf{W_i}\,(\mathrm{erg/s.cm}^{3})$',
						'figure/2D_maps/Wi/Wi_',
						0)

print('Background electron temperature Te')
dir= os.path.dirname("figure/2D_maps/Te/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/Te[eV].dat',
						Nx,Nz,
						'jet',
						r'$\log_{10}{\left(\mathbf{T_e}\,(\mathrm{eV})\right)}$',
						'figure/2D_maps/Te/Te_',
						1)

print('Background ion temperature Ti')
dir= os.path.dirname("figure/2D_maps/Ti/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/Ti[eV].dat',
						Nx,Nz,
						'jet',
						r'$\log_{10}{\left(\mathbf{T_i}\,(\mathrm{eV})\right)}$',
						'figure/2D_maps/Ti/Ti_',
						1)
						
print('Electrical resistivity eta')
dir= os.path.dirname("figure/2D_maps/eta/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/resis[Ohm.m].dat',
						Nx,Nz,
						'terrain',
						r'$\log_{10}{\left (\mathbf{\eta}\,(\Omega.\mathrm{m}) \right )}$',
						'figure/2D_maps/eta/eta_',
						1)
						
print('Thermal conductivity kappa')
dir= os.path.dirname("figure/2D_maps/kappa/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/Kappa_e[erg_m-1_K-1_s-1].dat',
						Nx,Nz,
						'terrain',
						r'$\log_{10}{\left (\mathbf{\kappa}\,(\mathrm{erg/m.K.s}) \right )}$',
						'figure/2D_maps/kappa/kappa_',
						1)

print('Ionization rate of K-shell electrons nuK')
dir= os.path.dirname("figure/2D_maps/nuK/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/K_shell_ioniz_rate_[s-1].dat',
						Nx,Nz,
						'OrRd',
						r'$\mathbf{\nu_K}\,(\mathrm{s}^{-1})$',
						'figure/2D_maps/nuK/nuK_',
						0)

print('Time integrated density of emitted Kalpha photons nKalpha')
dir= os.path.dirname("figure/2D_maps/nKa/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/n_Kalpha[cm-3].dat',
						Nx,Nz,
						'hot',
						r'$\mathbf{n_{K\alpha}}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1})$',
						'figure/2D_maps/nKa/nKa_',
						0)
					
print('Time integrated density of emitted Kbeta photons nKb')
dir= os.path.dirname("figure/2D_maps/nKb/")
if not os.path.exists(dir):
    os.mkdir(dir)
rapc.read_and_plot_2D_pcolormesh('diag/n_Kbeta[cm-3].dat',
						Nx,Nz,
						'hot',
						r'$\mathbf{n_{K\beta}}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1})$',
						'figure/2D_maps/nKb/nKb_',
						0)