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
Read and plot data X [X] | dN/dX [/[X]) from files:
* fast_electron_temporal_distr.dat : X = z [microns]
* fast_electron_spectrum.dat       : X = E [keV]
* fast_electron_spatial_distr.dat  : X = x [microns]
* fast_electron_angular_distr.dat  : X = theta [/degrees] but dN/dX [/rad]
"""
import numpy as np
import library as lib

SIMU_NAME=lib.get_results_dir()

print(' -------------------------------------')
print(' Electron Beam Initial Properties Plot')
print(' -------------------------------------')
print('  ')

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

INI_DIR = 'figures/'+SIMU_NAME+'/initialization/'
lib.create_dir(INI_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

print(' Injected electron beam kinetic energy spectrum')
print('  ')
RES_FILE  = RES_DIR+'electron_beam_kinetic_energy_spectrum.dat'
[E,DN_DE] = lib.read_file_and_define_two_first_cols(RES_FILE)
E_MIN     = lib.get_log_axis_min_value(E)
E_MAX     = lib.get_log_axis_max_value(E)
DN_DE_MIN = lib.get_log_axis_min_value(DN_DE)
DN_DE_MAX = lib.get_log_axis_max_value(DN_DE)
E_LBL     = r'$\varepsilon\,(\mathrm{keV})$'
DN_DE_LBL = r'$f_\varepsilon (\varepsilon)\,(/\mathrm{keV})$'
DN_DE_TTL = 'Injected electron beam kinetic energy spectrum'
DN_DE_FIG = INI_DIR+'electron_beam_kinetic_energy_spectrum.png'
lib.make_scalar_plot_figure(xplot     = E,
                            yplot     = DN_DE,
                            color     = 'red',
                            xplot_min = E_MIN,
                            xplot_max = E_MAX,
                            yplot_min = DN_DE_MIN,
                            yplot_max = DN_DE_MAX,
                            xlabel    = E_LBL,
                            ylabel    = DN_DE_LBL,
                            title     = DN_DE_TTL,
                            filename  = DN_DE_FIG,
                            logx      = True,
                            logy      = True,
                            grid      = True)

print(' Injected electron beam longitudinal distribution')
print('  ')
RES_FILE   = RES_DIR+'electron_beam_longitudinal_distribution.dat'
[T,DN_DZ]  = lib.read_file_and_define_two_first_cols(RES_FILE)
T_MAX      = np.max(T)
UNIT       = lib.get_log_axis_max_value(DN_DZ)
DN_DZ      = DN_DZ / UNIT
DN_DZ_MAX  = 1.1*np.max(DN_DZ)
T_LBL      = r'$t\,(\mathrm{fs})$'
DN_DZ_LBL  = r'$f_z \left (z=0,\,t \right )\,($'
DN_DZ_LBL += f"{UNIT:.0E}"+r'$/\mu\mathrm{m})$'
DN_DZ_TTL  = 'Injected electron beam longitudinal distribution'
DN_DZ_FIG  = INI_DIR+'electron_beam_longitudinal_distribution.png'
lib.make_scalar_plot_figure(xplot     = T,
                            yplot     = DN_DZ,
                            color     = 'red',
                            xplot_min = 0.,
                            xplot_max = T_MAX,
                            yplot_min = 0.,
                            yplot_max = DN_DZ_MAX,
                            xlabel    = T_LBL,
                            ylabel    = DN_DZ_LBL,
                            title     = DN_DZ_TTL,
                            filename  = DN_DZ_FIG,
                            logx      = False,
                            logy      = False,
                            grid      = False)

print(' Injected electron beam transverse distribution')
print('  ')
RES_FILE   = RES_DIR+'electron_beam_transverse_distribution.dat'
[X,DN_DX]  = lib.read_file_and_define_two_first_cols(RES_FILE)
X_MIN      = np.min(X)
X_MAX      = np.max(X)
UNIT       = lib.get_log_axis_max_value(DN_DX)
DN_DX      = DN_DX / UNIT
DN_DX_MAX  = 1.1*np.max(DN_DX)
X_LBL      = r'$x\,(\mu\mathrm{m})$'
DN_DX_LBL  = r'$f_x (x)\,($'
DN_DX_LBL += f"{UNIT:.0E}"+r'$/\mu\mathrm{m})$'
DN_DX_TTL  = 'Injected electron beam transverse distribution'
DN_DX_FIG  = INI_DIR+'electron_beam_transverse_distribution.png'
lib.make_scalar_plot_figure(xplot     = X,
                            yplot     = DN_DX,
                            color     = 'red',
                            xplot_min = X_MIN,
                            xplot_max = X_MAX,
                            yplot_min = 0.,
                            yplot_max = DN_DX_MAX,
                            xlabel    = X_LBL,
                            ylabel    = DN_DX_LBL,
                            title     = DN_DX_TTL,
                            filename  = DN_DX_FIG,
                            logx      = False,
                            logy      = False,
                            grid      = False)

print(' Injected electron beam angular distribution')
print('  ')

THETA_PLT = []
X_PLT     = []
P_PLT     = []
with open(RES_DIR+'electron_beam_angular_distribution.dat', 'r', encoding='utf-8') as file:
    for line in file:
        line      = line.strip()
        array     = line.split()
        X_PLT.append(float(array[0]))
        THETA_PLT.append(float(array[1]))
        P_PLT.append(float(array[2]))
N_THETA_N_X = len(X_PLT)
N_X         = len(X)
N_THETA     = int(N_THETA_N_X / N_X)
THETA_MAP   = np.zeros((N_THETA,N_X))
X_MAP       = np.zeros((N_THETA,N_X))
P_MAP       = np.zeros((N_THETA,N_X))
UNIT        = lib.get_log_axis_max_value(P_PLT)
P_PLT       = np.array(P_PLT) / UNIT
for i in range(0,N_X):
    for k in range(0,N_THETA):
        THETA_MAP[k][i] = THETA_PLT[i*N_THETA+k]
        X_MAP[k][i]     = X_PLT[i*N_THETA+k]
        P_MAP[k][i]     = P_PLT[i*N_THETA+k]
THETA_MIN      = np.min(THETA_PLT)
THETA_MAX      = np.max(THETA_PLT)
X_MIN          = np.min(X_PLT)
X_MAX          = np.max(X_PLT)
P_MIN          = np.min(P_PLT)
P_MAX          = np.max(P_PLT)
THETA_LBL      = r'$\theta\,(^\mathrm{o})$'
X_LBL          = r'$x\,(\mu\mathrm{m})$'
DN_DTHETA_TTL  = 'Injected electron beam angular distribution \n'
DN_DTHETA_TTL += r'$f_\theta \left (x,\,\theta \right )\,($'
DN_DTHETA_TTL += f"{UNIT:.0E}" + r'$/\mathrm{rad})$'
DN_DTHETA_FIG  = INI_DIR+'electron_beam_angular_distribution.png'
lib.make_2d_field_pcolormesh_figure(xmap     = THETA_MAP,
                                    ymap     = X_MAP,
                                    zmap     = P_MAP,
                                    colormap = 'jet',
                                    xmap_min = THETA_MIN,
                                    xmap_max = THETA_MAX,
                                    ymap_min = X_MIN,
                                    ymap_max = X_MAX,
                                    zmap_min = np.min(P_PLT),
                                    zmap_max = np.max(P_PLT),
                                    xlabel   = THETA_LBL,
                                    ylabel   = X_LBL,
                                    title    = DN_DTHETA_TTL,
                                    eq_axis  = False,
                                    fig_file = DN_DTHETA_FIG)
