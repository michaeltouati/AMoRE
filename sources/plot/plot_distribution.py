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
import numpy as np

SIMU_NAME=lib.get_results_dir()

print(' ------------------------------')
print(' Electron Beam Phase-space Plot')
print(' ------------------------------')

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

DIST_DIR = SIMU_DIR+'distribution_function/'
lib.create_dir(DIST_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

[N_Z,N_X]   = lib.find_spatial_simulation_box_dimension(RES_DIR+'B_y[Tesla].dat')

[N_EPS,EPS] = lib.find_energy_bins(RES_DIR+'psi0_x[cm-3_keV-1].dat')

N_THETA     = 500
THETA       = np.linspace(0., 2. * np.pi, N_THETA)
[DET_JACOBIAN, PZ_MAP, PX_MAP] = lib.make_phase_space_grid(n_theta = N_THETA,
                                                           theta   = THETA,
                                                           n_e     = N_EPS,
                                                           e_0     = EPS)
print('  ')
print(' At the location where the beam density is maximum at a given location on x-axis')
FBX_DIR = DIST_DIR+'fb_x/'
lib.create_dir(FBX_DIR)
PSI0_RES  = RES_DIR+'psi0_x[cm-3_keV-1].dat'
PSI1X_RES = RES_DIR+'psi1x_x[cm-3_keV-1].dat'
PSI1Z_RES = RES_DIR+'psi1z_x[cm-3_keV-1].dat'
FBX_FIG   = FBX_DIR+'fb_x'
lib.read_and_plot_distribution(n_theta        = N_THETA,
                               theta          = THETA,
                               n_e            = N_EPS,
                               e_0            = EPS,
                               det_jacobian   = DET_JACOBIAN,
                               pz_map         = PZ_MAP,
                               px_map         = PX_MAP,
                               mu_grk         = 2,
                               n_mu_grk       = N_X,
                               filename_psi0  = PSI0_RES,
                               filename_psi1x = PSI1X_RES,
                               filename_psi1z = PSI1Z_RES,
                               fig_name       = FBX_FIG)
print('  ')
print(' At the location where the beam density is maximum at a given depth')
FBZ_DIR = DIST_DIR+'fb_z/'
lib.create_dir(FBZ_DIR)
PSI0_RES  = RES_DIR+'psi0_z[cm-3_keV-1].dat'
PSI1X_RES = RES_DIR+'psi1x_z[cm-3_keV-1].dat'
PSI1Z_RES = RES_DIR+'psi1z_z[cm-3_keV-1].dat'
FBZ_FIG   = FBZ_DIR+'fb_z'
lib.read_and_plot_distribution(n_theta        = N_THETA,
                               theta          = THETA,
                               n_e            = N_EPS,
                               e_0            = EPS,
                               det_jacobian   = DET_JACOBIAN,
                               pz_map         = PZ_MAP,
                               px_map         = PX_MAP,
                               mu_grk         = 1,
                               n_mu_grk       = N_Z,
                               filename_psi0  = PSI0_RES,
                               filename_psi1x = PSI1X_RES,
                               filename_psi1z = PSI1Z_RES,
                               fig_name       = FBZ_FIG)
