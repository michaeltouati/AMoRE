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
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import rapc

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)
    							
t       = []
rapc.read_file_and_define_first_col('diag/U_e[J].dat',t)
Ue      = []
rapc.read_file_and_define_second_col('diag/U_e[J].dat',Ue)
Ub      = []
rapc.read_file_and_define_second_col('diag/U_b[J].dat',Ub)
Ud_col  = []
rapc.read_file_and_define_second_col('diag/Ud_col[J].dat', Ud_col)
Ud_res  = []
rapc.read_file_and_define_second_col('diag/Ud_res[J].dat', Ud_res)
U_el    = []
rapc.read_file_and_define_second_col('diag/U_el[J].dat',   U_el)
U_ma    = []
rapc.read_file_and_define_second_col('diag/U_ma[J].dat',   U_ma)
U_sd    = []
rapc.read_file_and_define_second_col('diag/U_sd[J].dat',   U_sd)
U_su    = []
rapc.read_file_and_define_second_col('diag/U_su[J].dat',   U_su)
U_sf    = []
rapc.read_file_and_define_second_col('diag/U_sf[J].dat',   U_sf)
U_sb    = []
rapc.read_file_and_define_second_col('diag/U_sb[J].dat',   U_sb)
if np.amax(U_el)>0:
	coef_el = np.floor(np.log(np.amax(Ub)/np.amax(U_el))/np.log(10))
	coef_el = 10.**coef_el
else:
	coef_el = 1.
coef_ma = np.floor(np.log(np.amax(Ub)/np.amax(U_ma))/np.log(10))
coef_ma = 10.**coef_ma
U_stot  = []
U_el_n  = []
U_ma_n  = []
nn=len(t)-1
for i in range(0, nn+1):
	U_stot.append(U_su[i]+U_sd[i]+U_sf[i]+U_sb[i])
	U_el_n.append(coef_el*U_el[i])
	U_ma_n.append(coef_ma*U_ma[i])
error = 100.*( Ue[nn]-(Ub[nn]+U_el[nn]+U_ma[nn]+Ud_col[nn]+Ud_res[nn]+U_stot[nn]) )/Ue[nn]
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(t, U_el_n,'green',linewidth=2,label=r'$\mathbf{U_{E}\times}$'+str('%.0e' % coef_el))
plt.plot(t, U_ma_n,'cyan',linewidth=2,label=r'$\mathbf{U_{B}\times}$'+str('%.0e' % coef_ma))
plt.plot(t, Ue,'red',linewidth=2,label=r'$\mathbf{U_\mathrm{inj}}$')
plt.plot(t, Ub,'black',linewidth=2,label=r'$\mathbf{U_b}$')
plt.plot(t, Ud_col,'magenta',linewidth=2,label=r'$\mathbf{U_{d,\mathrm{col}}}$')
plt.plot(t, Ud_res,'blue',linewidth=2,label=r'$\mathbf{U_{d,\mathrm{res}}}$')
plt.plot(t, U_stot,'red',linestyle='--',linewidth=2,label=r'$\mathbf{U_\mathrm{esc}}$')
leg = plt.legend(loc='upper left',fontsize=16, fancybox=True)
leg.get_frame().set_alpha(0.5)
#font = {'family': 'non-serif',
font = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.title('Energy conservation error '+'='+str(np.floor(100*error)/100)+r'$\,\%$', fontdict=font)
plt.xticks(fontsize=16)
plt.xlabel('time (fs)', fontdict=font)
plt.xlim([t[0],t[nn]])
plt.ylabel('Energy (J)', fontdict=font)
plt.yticks(fontsize=16)
plt.ylim([Ue[0],1.1*Ue[nn]])
fig.savefig('figures/energy_conservation.png',bbox_inches='tight')
plt.close(fig)
