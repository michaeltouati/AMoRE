import numpy as np
import math
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import rapc

dir= os.path.dirname("figure/")
if not os.path.exists(dir):
    os.mkdir(dir)

t     = []
rapc.read_file_and_define_first_col('diag/fast_electron_temporal_distr.dat',t)
dN_dt   = []
rapc.read_file_and_define_second_col('diag/fast_electron_temporal_distr.dat',dN_dt)

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
fig.savefig('figure/injected_fast_e_temporal_distr.png',bbox_inches='tight')
plt.close(fig)



eps     = []
rapc.read_file_and_define_first_col('diag/fast_electron_spectrum.dat',eps)
dN_dE   = []
rapc.read_file_and_define_second_col('diag/fast_electron_spectrum.dat',dN_dE)

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
fig.savefig('figure/injected_fast_e_energy_spectrum.png',bbox_inches='tight')
plt.close(fig)

x     = []
rapc.read_file_and_define_first_col('diag/fast_electron_spatial_distr.dat',x)
dN_dx   = []
rapc.read_file_and_define_second_col('diag/fast_electron_spatial_distr.dat',dN_dx)

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
fig.savefig('figure/injected_fast_e_transverse_distribution.png',bbox_inches='tight')
plt.close(fig)

N1 = 360
N2 = len(x)
N3 = N1*N2
Theta = np.zeros((N1,N2))
X     = np.zeros((N1,N2))
P     = np.zeros((N1,N2))
file = open('diag/fast_electron_angular_distr.dat', 'r')
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
fig.savefig('figure/injected_fast_e_angular_distribution.png',bbox_inches='tight')
plt.close(fig)

