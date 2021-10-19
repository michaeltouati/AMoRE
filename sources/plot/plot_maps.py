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
import library as lib
import os

simu_name=lib.get_results_dir()

print(' ---------------------------------------')
print(' Hydrodynamic Moments 2D Density Plot')
print(' ---------------------------------------')
print('  ')

lib.create_dir('figures/')
lib.create_dir('figures/'+simu_name+'/')

subdir = 'figures/'+simu_name+'/2D_maps/'
lib.create_dir(subdir)

results_dir = 'results/'+simu_name+'/'

[N_z,N_x] = lib.find_spatial_simulation_box_dimension(results_dir+'B_y[Tesla].dat')

# plot 2D maps

print(' Target electron density ne')
subdir2 = subdir+'ne/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'ne[cm-3].dat',
								N_x,N_z,
								'jet',
								r'$\mathbf{n_e}\,(\mathrm{cm}^{-3})$',
								subdir2+'ne_',
								0)
print('  ')
print(' Beam density nb')
subdir2 = subdir+'nb/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'nb[cm-3].dat',
								N_x,N_z,
								'Blues',
								r'$\mathbf{n_b}\,(\mathrm{cm}^{-3})$',
								subdir2+'nb_',
								0)
print('  ')
print(' Target ion density ni')
subdir2 = subdir+'ni/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'ni[cm-3].dat',
								N_x,N_z,
								'jet',
								r'$\mathbf{n_i}\,(\mathrm{cm}^{-3})$',
								subdir2+'ni_',
								0)
print('  ')
print(' Beam current density jbz')
subdir2 = subdir+'jbz/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'jb_z[A_cm-2].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{j_{b,z}}\,(\mathrm{A/cm}^{2})$',
								subdir2+'jbz_',
								0)
# print('  ')
print(' Beam current density jbx') 
subdir2 = subdir+'jbx/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'jb_x[A_cm-2].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{j_{b,x}}\,(\mathrm{A/cm}^{2})$',
								subdir2+'jbx_',
								0)
# print('  ')
# print(' Beam current density |jb|')    
# subdir2 = subdir+'jb/'
# lib.create_dir(subdir2)
# lib.read_and_plot_2D_pcolormesh_abs(results_dir+'jb_z[A_cm-2].dat',
# 									results_dir+'jb_x[A_cm-2].dat',
# 									N_x,N_z,
# 									'Reds',
# 									r'$|\mathbf{j_b}\,(\mathrm{A/cm}^{2})|$',
# 									subdir2+'jb_')
print('  ')
print(' Return current density jez')
subdir2 = subdir+'jez/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'je_z[A_cm-2].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{j_{e,z}}\,(\mathrm{A/cm}^{2})$',
								subdir2+'jez_',
								0)
print('  ')
print(' Return current density jex') 
subdir2 = subdir+'jex/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'je_x[A_cm-2].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{j_{e,x}}\,(\mathrm{A/cm}^{2})$',
								subdir2+'jex_',
								0)
# print('  ')
# print(' Return current density |je|')    
# subdir2 = subdir+'je/'
# lib.create_dir(subdir2)
# lib.read_and_plot_2D_pcolormesh_abs(results_dir+'je_z[A_cm-2].dat',
# 									results_dir+'je_x[A_cm-2].dat',
# 									N_x,N_z,
# 									'Reds',
# 									r'$|\mathbf{j_e}\,(\mathrm{A/cm}^{2})|$',
# 									subdir2+'je_')
print('  ')
print(' Electric field Ex')
subdir2 = subdir+'Ex/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'E_x[V_m-1].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{E_{x}}\,(\mathrm{V/m})$',
								subdir2+'Ex_',
								0)
print('  ')
print(' Electric field Ez')
subdir2 = subdir+'Ez/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'E_z[V_m-1].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{E_{z}}\,(\mathrm{V/m})$',
								subdir2+'Ez_',
								0)
# print('  ')
# print(' Electric field |E|')    
# subdir2 = subdir+'E/'
# lib.create_dir(subdir2)
# lib.read_and_plot_2D_pcolormesh_abs(results_dir+'E_z[V_m-1].dat',
# 									results_dir+'E_x[V_m-1].dat',
# 									N_x,N_z,
# 									'Reds',
# 									r'$|\mathbf{E}\,(\mathrm{V/m})|$',
# 									subdir2+'E_')
print('  ')
print(' Magnetic field By')
subdir2 = subdir+'By/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'B_y[Tesla].dat',
								N_x,N_z,
								'seismic',
								r'$\mathbf{B_y}\,(\mathrm{T})$',
								subdir2+'By_',
								0)
print('  ')
print(' Density of power deposited on background electrons We')
subdir2 = subdir+'We/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'We[erg_s-1_cm-3].dat',
								N_x,N_z,
								'OrRd',
								r'$\mathbf{W_e}\,(\mathrm{erg/s.cm}^{3})$',
								subdir2+'We_',
								0)
print('  ')						
print(' Density of power deposited on background ions Wi')
subdir2 = subdir+'Wi/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'Wi[erg_s-1_cm-3].dat',
								N_x,N_z,
								'OrRd',
								r'$\mathbf{W_i}\,(\mathrm{erg/s.cm}^{3})$',
								subdir2+'Wi_',
								0)
print('  ')
print(' Background electron temperature Te')
subdir2 = subdir+'Te/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'Te[eV].dat',
								N_x,N_z,
								'jet',
								r'$\log_{10}{\left(\mathbf{T_e}\,(\mathrm{eV})\right)}$',
								subdir2+'Te_',
								1)
print('  ')
print(' Background ion temperature Ti')
subdir2 = subdir+'Ti/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'Ti[eV].dat',
								N_x,N_z,
								'jet',
								r'$\log_{10}{\left(\mathbf{T_i}\,(\mathrm{eV})\right)}$',
								subdir2+'Ti_',
								1)
print('  ')						
print(' Electrical resistivity eta')
subdir2 = subdir+'eta/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'resis[Ohm.m].dat',
								N_x,N_z,
								'terrain',
								r'$\log_{10}{\left (\mathbf{\eta}\,(\Omega.\mathrm{m}) \right )}$',
								subdir2+'eta_',
								1)
print('  ')						
print(' Thermal conductivity kappa')
subdir2 = subdir+'kappa/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'Kappa_e[erg_m-1_K-1_s-1].dat',
								N_x,N_z,
								'terrain',
								r'$\log_{10}{\left (\mathbf{\kappa}\,(\mathrm{erg/m.K.s}) \right )}$',
								subdir2+'kappa_',
								1)
print('  ')
print(' Ionization rate of K-shell electrons nuK')
subdir2 = subdir+'nuK/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'K_shell_ioniz_rate_[s-1].dat',
								N_x,N_z,
								'OrRd',
								r'$\mathbf{\nu_K}\,(\mathrm{s}^{-1})$',
								subdir2+'nuK_',
								0)
print('  ')
print(' Time integrated density of emitted Kalpha photons nKalpha')
subdir2 = subdir+'nKa/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'n_Kalpha[cm-3].dat',
								N_x,N_z,
								'hot',
								r'$\mathbf{n_{K\alpha}}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1})$',
								subdir2+'nKa_',
								0)
print('  ')
print(' Time integrated density of emitted Kbeta photons nKb')
subdir2 = subdir+'nKb/'
lib.create_dir(subdir2)
lib.read_and_plot_2D_pcolormesh(results_dir+'n_Kbeta[cm-3].dat',
								N_x,N_z,
								'hot',
								r'$\mathbf{n_{K\beta}}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1})$',
								subdir2+'nKb_',
								0)
print('  ')