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
"""
Read and plot the electron beam distribution function from its
two first momentum angular moments files t (fs) | x (microns) | z (microns) | F :
* psi0_x[cm-3_keV-1].dat
* psi1x_x[cm-3_keV-1].dat or psi1x_z[cm-3_keV-1].dat
* psi1z_x[cm-3_keV-1].dat or psi1z_z[cm-3_keV-1].dat
"""
import library as lib

simu_name=lib.get_results_dir()

print(' ------------------------------')
print(' Electron Beam Phase-space Plot')
print(' ------------------------------')
print('  ')

lib.create_dir('figures/')
lib.create_dir('figures/'+simu_name+'/')

subdir = 'figures/'+simu_name+'/distribution_function/'
lib.create_dir(subdir)

results_dir = 'results/'+simu_name+'/'

[N_z,N_x] = lib.find_spatial_simulation_box_dimension(results_dir+'B_y[Tesla].dat')

[N_eps,eps] = lib.find_energy_bins(results_dir+'psi0_x[cm-3_keV-1].dat')

print(' At the location where the beam density is maximum at a given location on x-axis')
print('  ')
subdir2 = subdir+'fb_x/'
lib.create_dir(subdir2)
lib.read_and_plot_distribution(n_e=N_eps,
                               e_0=eps,
                               mu_grk=2,
                               n_mu_grk=N_x,
                               filename_psi0=results_dir+'psi0_x[cm-3_keV-1].dat',
                               filename_psi1x=results_dir+'psi1x_x[cm-3_keV-1].dat',
                               filename_psi1z=results_dir+'psi1z_x[cm-3_keV-1].dat',
                               name=subdir2+'fb_x')

print(' At the location where the beam density is maximum at a given depth')
print('  ')
subdir2 = subdir+'fb_z/'
lib.create_dir(subdir2)
lib.read_and_plot_distribution(n_e=N_eps,
                               e_0=eps,
                               mu_grk=1,
                               n_mu_grk=N_z,
                               filename_psi0=results_dir+'psi0_z[cm-3_keV-1].dat',
                               filename_psi1x=results_dir+'psi1x_z[cm-3_keV-1].dat',
                               filename_psi1z=results_dir+'psi1z_z[cm-3_keV-1].dat',
                               name=subdir2+'fb_z')
