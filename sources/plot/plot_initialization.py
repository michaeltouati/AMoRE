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

print(' Injected Electron Beam Longitudinal Distribution')
print('  ')
RES_FILE  = RES_DIR+'electron_beam_longitudinal_distribution.dat'
[T,DN_DZ] = lib.read_file_and_define_two_first_cols(RES_FILE)
T_MAX     = np.max(T)
UNIT      = lib.get_log_axis_max_value(DN_DZ)
DN_DZ     = DN_DZ / UNIT
DN_DZ_MAX = 1.1*np.max(DN_DZ)
Y_LABEL   = r'$f_z \left (z=0,\,t \right )\,($'
Y_LABEL   = Y_LABEL +f"{UNIT:.0E}"+r'$/\mu\mathrm{m})$'
FIG_FILE  = INI_DIR+'electron_beam_longitudinal_distribution.png'
lib.make_scalar_plot_figure(xplot=T, yplot=DN_DZ,
                            color='red',
                            xplot_min=0.,
                            xplot_max=T_MAX,
                            yplot_min=0.,
                            yplot_max=DN_DZ_MAX,
                            xlabel=r'$t\,(\mathrm{fs})$',
                            ylabel=Y_LABEL,
                            title='Injected Electron Beam Longitudinal Distribution',
                            filename=FIG_FILE,
                            logx=False, logy=False, grid=False)

print(' Injected Electron Beam Kinetic Energy Spectrum')
print('  ')
RES_FILE  = RES_DIR+'electron_beam_kinetic_energy_spectrum.dat'
[E,DN_DE] = lib.read_file_and_define_two_first_cols(RES_FILE)
E_MIN     = lib.get_log_axis_min_value(E)
E_MAX     = lib.get_log_axis_max_value(E)
DN_DE_MIN = lib.get_log_axis_min_value(DN_DE)
DN_DE_MAX = lib.get_log_axis_max_value(DN_DE)
FIG_FILE  = INI_DIR+'electron_beam_kinetic_energy_spectrum.png'
lib.make_scalar_plot_figure(xplot=E, yplot=DN_DE,
                            color='red',
                            xplot_min=E_MIN,
                            xplot_max=E_MAX,
                            yplot_min=DN_DE_MIN,
                            yplot_max=DN_DE_MAX,
                            xlabel=r'$\varepsilon\,(\mathrm{keV})$',
                            ylabel=r'$f_\varepsilon (\varepsilon)\,(/\mathrm{keV})$',
                            title='Injected Electron Beam Kinetic Energy Spectrum',
                            filename=FIG_FILE,
                            logx=True, logy=True, grid=True)

print(' Injected Electron Beam Transverse Distribution')
print('  ')
RES_FILE  = RES_DIR+'electron_beam_transverse_distribution.dat'
[X,DN_DX] = lib.read_file_and_define_two_first_cols(RES_FILE)
X_MIN     = np.min(X)
X_MAX     = np.max(X)
UNIT      = lib.get_log_axis_max_value(DN_DX)
DN_DX     = DN_DX / UNIT
DN_DX_MAX = 1.1*np.max(DN_DX)
Y_LABEL = r'$f_x (x)\,($'
Y_LABEL = Y_LABEL +f"{UNIT:.0E}"+r'$/\mu\mathrm{m})$'
FIG_FILE  = INI_DIR+'electron_beam_transverse_distribution.png'
lib.make_scalar_plot_figure(xplot=X, yplot=DN_DX,
                            color='red',
                            xplot_min=X_MIN,
                            xplot_max=X_MAX,
                            yplot_min=0.,
                            yplot_max=DN_DX_MAX,
                            xlabel=r'$x\,(\mu\mathrm{m})$',
                            ylabel=Y_LABEL,
                            title='Injected Electron Beam Transverse Distribution',
                            filename=FIG_FILE,
                            logx=False, logy=False, grid=False)

print(' Injected Electron Beam Angular Distribution')
print('  ')

N1 = 360
N2 = len(X)
N3 = N1*N2
Theta = np.zeros((N1,N2))
X     = np.zeros((N1,N2))
P     = np.zeros((N1,N2))
theta = []
x2    = []
p     = []
with open(RES_DIR+'electron_beam_angular_distribution.dat', 'r', encoding='utf-8') as file:
    for line in file:
        line      = line.strip()
        array     = line.split()
        x2.append(float(array[0]))
        theta.append(float(array[1]))
        p.append(float(array[2]))
UNIT   = lib.get_log_axis_max_value(p)
p      = np.array(p) / UNIT
Maxval = np.max(p)
Minval = np.min(p)
for i in range(0,N2):
    for k in range(0,N1):
        Theta[k][i]=theta[i*N1+k]
        X[k][i]    =x2[i*N1+k]
        P[k][i]    =p[i*N1+k]
TTL  = 'Injected Electron Beam Angular Distribution '
TTL += r'$f_\theta \left (x,\,\theta \right )\,($'
TTL += f"{UNIT:.0E}" + r'$/\mathrm{rad})$'
FIG_FILE = INI_DIR+'electron_beam_angular_distribution.png'
lib.make_2d_field_pcolormesh_figure(xmap=Theta,
                                    ymap=X,
                                    zmap=P,
                                    colormap='jet',
                                    xmap_min=np.min(theta),
                                    xmap_max=np.max(theta),
                                    ymap_min=np.min(X),
                                    ymap_max=np.max(X),
                                    zmap_min=np.min(p),
                                    zmap_max=np.max(p),
                                    xlabel=r'$\theta\,(^\mathrm{o})$',
                                    ylabel=r'$x\,(\mu\mathrm{m})$',
                                    title=TTL,
                                    filename=FIG_FILE)
