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
Read and plot data t (fs) | Energy (J) from files:
* U_b[J].dat
* U_e[J].dat
* U_el[J].dat
* U_ma[J].dat
* U_sb[J].dat
* U_sd[J].dat
* U_sf[J].dat
* U_su[J].dat
"""
import matplotlib.pyplot as plt
import numpy as np
import library as lib

simu_name=lib.get_results_dir()

print(' --------------------')
print(' Energies Scalar Plot')
print(' --------------------')
print('  ')

lib.create_dir('figures/')
subdir = 'figures/'+simu_name+'/'
lib.create_dir(subdir)

results_dir = 'results/'+simu_name+'/'

T      = lib.read_file_and_define_first_col(results_dir+'U_e[J].dat')
U_E    = lib.read_file_and_define_second_col(results_dir+'U_e[J].dat')
U_B    = lib.read_file_and_define_second_col(results_dir+'U_b[J].dat')
U_Dcol = lib.read_file_and_define_second_col(results_dir+'Ud_col[J].dat')
U_Dres = lib.read_file_and_define_second_col(results_dir+'Ud_res[J].dat')
U_EL   = lib.read_file_and_define_second_col(results_dir+'U_el[J].dat')
U_MA   = lib.read_file_and_define_second_col(results_dir+'U_ma[J].dat')
U_Sd   = lib.read_file_and_define_second_col(results_dir+'U_sd[J].dat')
U_Su   = lib.read_file_and_define_second_col(results_dir+'U_su[J].dat')
U_Sf   = lib.read_file_and_define_second_col(results_dir+'U_sf[J].dat')
U_Sb   = lib.read_file_and_define_second_col(results_dir+'U_sb[J].dat')
if np.amax(U_EL)>0:
    C_EL = np.floor(np.log(np.amax(U_B)/np.amax(U_EL))/np.log(10))
    C_EL = 10.**C_EL
else:
    C_EL = 1.
C_MA = np.floor(np.log(np.amax(U_B)/np.amax(U_MA))/np.log(10))
C_MA = 10.**C_MA
U_stot  = []
U_EL_n  = []
U_MA_n  = []
N_T=len(T)-1
for i in range(0, N_T+1):
    U_stot.append(U_Su[i]+U_Sd[i]+U_Sf[i]+U_Sb[i])
    U_EL_n.append(C_EL*U_EL[i])
    U_MA_n.append(C_MA*U_MA[i])
balance = U_E[N_T]-(U_B[N_T]+U_EL[N_T]+U_MA[N_T]+U_Dcol[N_T]+U_Dres[N_T]+U_stot[N_T])
error   = 100.*balance/U_E[N_T]
fig=plt.figure()
plt.rc('text', usetex=True)
plt.plot(T, U_EL_n,'green',linewidth=2,label=r'$\mathbf{U_{E}\times}$'+f"{C_EL:.0E}")
plt.plot(T, U_MA_n,'cyan',linewidth=2,label=r'$\mathbf{U_{B}\times}$'+f"{C_MA:.0E}")
plt.plot(T, U_E,'red',linewidth=2,label=r'$\mathbf{U_\mathrm{inj}}$')
plt.plot(T, U_B,'black',linewidth=2,label=r'$\mathbf{U_b}$')
plt.plot(T, U_Dcol,'magenta',linewidth=2,label=r'$\mathbf{U_{d,\mathrm{col}}}$')
plt.plot(T, U_Dres,'blue',linewidth=2,label=r'$\mathbf{U_{d,\mathrm{res}}}$')
plt.plot(T, U_stot,'red',linestyle='--',linewidth=2,label=r'$\mathbf{U_\mathrm{esc}}$')
leg = plt.legend(loc='upper left',fontsize=lib.font_size, fancybox=True)
leg.get_frame().set_alpha(0.5)
TTL = 'Energy conservation error '
TTL = TTL + '='+str(np.floor(100*error)/100)+r'$\,\%$'
plt.title(TTL, fontdict=lib.font)
plt.xticks(fontsize=lib.font_size)
plt.xlabel('time (fs)', fontdict=lib.font)
plt.xlim([T[0],T[N_T]])
plt.ylabel('Energy (J)', fontdict=lib.font)
plt.yticks(fontsize=lib.font_size)
plt.ylim([U_E[0],1.1*U_E[N_T]])
fig.savefig(subdir+'energy_conservation.png',bbox_inches='tight')
plt.close(fig)
