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
import matplotlib.pyplot as plt
from matplotlib import cm

FONT_SIZE=16

FONT = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': FONT_SIZE,
        }

simu_name=lib.get_results_dir()

print(' -------------------------------------')
print(' Electron Beam Initial Properties Plot')
print(' -------------------------------------')
print('  ')

lib.create_dir('figures/')
lib.create_dir('figures/'+simu_name+'/')

subdir = 'figures/'+simu_name+'/initialization/'
lib.create_dir(subdir)

results_dir = 'results/'+simu_name+'/'

print(' Injected Electron Beam Longitudinal Distribution')
print('  ')

[t,dN_dz] = lib.read_file_and_define_two_first_cols(results_dir+'fast_electron_temporal_distr.dat')
maxval    = np.max(dN_dz)
unit      = 10.**(1.+np.int(np.log(maxval)/np.log(10.)))
dN_dz     = np.array(dN_dz) / unit
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(t, dN_dz,'black',linewidth=2)
TTL = 'Injected Fast Electron Longitudinal Distribution'
Y_LABEL = r'$f_z \left (z=0,\,t \right )\,($'
Y_LABEL = Y_LABEL +f"{unit:.0E}"+r'$/\mu\mathrm{m})$'
plt.title(TTL, fontdict=FONT)
plt.xticks(fontsize=FONT_SIZE)
plt.xlabel(r'$t\,(\mathrm{fs})$', fontdict=FONT)
plt.xlim([np.min(t),np.max(t)])
plt.ylabel(Y_LABEL, fontdict=FONT)
plt.yticks(fontsize=FONT_SIZE)
plt.ylim([np.min(dN_dz),1.1*np.max(dN_dz)])
fig.savefig(subdir+'injected_fast_e_temporal_distr.png',bbox_inches='tight')
plt.close(fig)

print(' Injected Electron Beam Kinetic Energy Spectrum')
print('  ')

[eps,dN_dE] = lib.read_file_and_define_two_first_cols(results_dir+'fast_electron_spectrum.dat')
maxval    = np.max(dN_dE)
unit      = 10.**(1.+np.int(np.log(maxval)/np.log(10.)))
dN_dE     = np.array(dN_dE) / unit
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(eps, dN_dE,'black',linewidth=2)
TTL = 'Injected Fast Electron Energy Spectrum'
plt.title(TTL, fontdict=FONT)
plt.xticks(fontsize=FONT_SIZE)
plt.xlabel(r'$\varepsilon\,(\mathrm{keV})$', fontdict=FONT)
plt.xlim([np.min(eps),np.max(eps)])
Y_LABEL = r'$f_\varepsilon (\varepsilon)\,($'
Y_LABEL = Y_LABEL +f"{unit:.0E}"+r'$/\mathrm{keV})$'
plt.ylabel(Y_LABEL, fontdict=FONT)
plt.yticks(fontsize=FONT_SIZE)
plt.ylim([np.min(dN_dE),1.1*np.max(dN_dE)])
fig.savefig(subdir+'injected_fast_e_energy_spectrum.png',bbox_inches='tight')
plt.close(fig)

print(' Injected Electron Beam Transverse Distribution')
print('  ')

[x,dN_dx] = lib.read_file_and_define_two_first_cols(results_dir+'fast_electron_spatial_distr.dat')
maxval    = np.max(dN_dx)
unit      = 10.**(1.+np.int(np.log(maxval)/np.log(10.)))
dN_dx     = np.array(dN_dx) / unit
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(x, dN_dx,'black',linewidth=2)
TTL = 'Injected Fast Electron Transverse Distribution'
plt.title(TTL, fontdict=FONT)
plt.xticks(fontsize=FONT_SIZE)
plt.xlabel(r'$x\,(\mu\mathrm{m})$', fontdict=FONT)
plt.xlim([np.min(x),np.max(x)])
Y_LABEL = r'$f_x (x)\,($'
Y_LABEL = Y_LABEL +f"{unit:.0E}"+r'$/\mu\mathrm{m})$'
plt.ylabel(Y_LABEL, fontdict=FONT)
plt.yticks(fontsize=FONT_SIZE)
plt.ylim([np.min(dN_dx),1.1*np.max(dN_dx)])
fig.savefig(subdir+'injected_fast_e_transverse_distribution.png',bbox_inches='tight')
plt.close(fig)

print(' Injected Electron Beam Angular Distribution')
print('  ')

N1 = 360
N2 = len(x)
N3 = N1*N2
Theta = np.zeros((N1,N2))
X     = np.zeros((N1,N2))
P     = np.zeros((N1,N2))
theta = []
x2    = []
p     = []
with open(results_dir+'fast_electron_angular_distr.dat', 'r', encoding='utf-8') as file:
    for line in file:
        line      = line.strip()
        array     = line.split()
        x2.append(float(array[0]))
        theta.append(float(array[1]))
        p.append(float(array[2]))
Maxval = np.max(p)
Minval = np.min(p)
unit   = 10.**(1.+np.int(np.log(maxval)/np.log(10.)))
p      = np.array(p) / unit
Maxval = Maxval / unit
Minval = Minval / unit
for i in range(0,N2):
    for k in range(0,N1):
        Theta[k][i]=theta[i*N1+k]
        X[k][i]    =x2[i*N1+k]
        P[k][i]    =p[i*N1+k]
cmap = plt.get_cmap('jet')
norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
fig=plt.figure()
plt.rc('text', usetex=True)
plt.pcolormesh(Theta,X,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
cbar=plt.colorbar()
cbar.ax.tick_params(labelsize=FONT_SIZE)
TTL = 'Injected Fast Electron Angular Distribution '
TTL = TTL + r'$f_\theta \left (x,\,\theta \right )\,($'
TTL = TTL + f"{unit:.0E}" + r'$/\mathrm{rad})$'
plt.title(TTL, fontdict=FONT)
plt.xticks(fontsize=16)
plt.xlabel(r'$\theta\,(^\mathrm{o})$', fontdict=FONT)
plt.xlim([np.min(theta),np.max(theta)])
plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=FONT)
plt.yticks(fontsize=16)
plt.ylim([np.min(x2),np.max(x2)])
fig.savefig(subdir+'injected_fast_e_angular_distribution.png',bbox_inches='tight')
plt.close(fig)
