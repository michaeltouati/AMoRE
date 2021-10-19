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
import os
import math
import numpy as np
import library as lib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

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

t     = []
lib.read_file_and_define_first_col(results_dir+'fast_electron_temporal_distr.dat',t)
dN_dt   = []
lib.read_file_and_define_second_col(results_dir+'fast_electron_temporal_distr.dat',dN_dt)

print(' Injected Electron Beam Longitudinal Distribution')

fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(t, dN_dt,'black',linewidth=2)
#font = {'family': 'non-serif',
font = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.title('Injected Fast Electron Temporal Distribution', fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel(r'$t\,(\mathrm{fs})$', fontdict=font)
plt.xlim([np.min(t),np.max(t)])
plt.ylabel(r'$d N / d t\,(\mathrm{a.u.})$', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([np.min(dN_dt),1.1*np.max(dN_dt)])
fig.savefig(subdir+'injected_fast_e_temporal_distr.png',bbox_inches='tight')
plt.close(fig)
print('  ')
print(' Injected Electron Beam Kinetic Energy Spectrum')

eps     = []
lib.read_file_and_define_first_col(results_dir+'fast_electron_spectrum.dat',eps)
dN_dE   = []
lib.read_file_and_define_second_col(results_dir+'fast_electron_spectrum.dat',dN_dE)

fig=plt.figure()
plt.rc('text', usetex=True)
#plt.loglog(eps, dN_dE,'black',linewidth=2)
plt.plot(eps, dN_dE,'black',linewidth=2)
#font = {'family': 'non-serif',
font = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.title('Injected Fast Electron Energy Spectrum', fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel(r'$\varepsilon\,(\mathrm{keV})$', fontdict=font)
plt.xlim([np.min(eps),np.max(eps)])
plt.ylabel(r'$d N / d \varepsilon\,(\mathrm{a.u.})$', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([np.min(dN_dE),1.1*np.max(dN_dE)])
fig.savefig(subdir+'injected_fast_e_energy_spectrum.png',bbox_inches='tight')
plt.close(fig)

x     = []
lib.read_file_and_define_first_col(results_dir+'fast_electron_spatial_distr.dat',x)
dN_dx   = []
lib.read_file_and_define_second_col(results_dir+'fast_electron_spatial_distr.dat',dN_dx)
print('  ')
print(' Injected Electron Beam Transverse Distribution')

fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(x, dN_dx,'black',linewidth=2)
#font = {'family': 'non-serif',
font = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.title('Injected Fast Electron Transverse Distribution', fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel(r'$x\,(\mu\mathrm{m})$', fontdict=font)
plt.xlim([np.min(x),np.max(x)])
plt.ylabel(r'$d N / d x\,(\mathrm{a.u.})$', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([np.min(dN_dx),1.1*np.max(dN_dx)])
fig.savefig(subdir+'injected_fast_e_transverse_distribution.png',bbox_inches='tight')
plt.close(fig)
print('  ')
print(' Injected Electron Beam Angular Distribution')

N1 = 360
N2 = len(x)
N3 = N1*N2
Theta = np.zeros((N1,N2))
X     = np.zeros((N1,N2))
P     = np.zeros((N1,N2))
file = open(results_dir+'fast_electron_angular_distr.dat', 'r')
theta = []
x2    = []
p     = [] 
for line in file:
	line      = line.strip()
	array     = line.split()
	x2.append(float(array[0]))
	theta.append(float(array[1]))
	p.append(float(array[2]))
file.close()
for i in range(0,N2):
	for k in range(0,N1):
		Theta[k][i]=theta[i*N1+k]
		X[k][i]    =x2[i*N1+k]
		P[k][i]    =p[i*N1+k]
cmap = plt.get_cmap('jet')
Maxval = np.max(p)
Minval = np.min(p)
norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
fig=plt.figure()
plt.rc('text', usetex=True)
plt.pcolormesh(Theta,X,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
cbar=plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.gca().set_aspect('equal')
#font = {'family': 'non-serif',
font = {'style':  'normal',
        'color':  'black',
       	'weight': 'normal',
       	'size': 16,
       	}
plt.title('Injected Fast Electron Angular Distribution', fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel(r'$\theta\,(^\mathrm{o})$', fontdict=font)
plt.xlim([np.min(theta),np.max(theta)])
plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([np.min(x2),np.max(x2)])
fig.savefig(subdir+'injected_fast_e_angular_distribution.png',bbox_inches='tight')
plt.close(fig)

print('  ')