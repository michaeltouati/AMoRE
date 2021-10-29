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
import numpy as np
import library as lib

SIMU_NAME=lib.get_results_dir()

print(' --------------------')
print(' Energies Scalar Plot')
print(' --------------------')
print('  ')

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

T      = lib.read_file_and_define_first_col( RES_DIR+'U_e[J].dat')
U_E    = lib.read_file_and_define_second_col(RES_DIR+'U_e[J].dat')
U_B    = lib.read_file_and_define_second_col(RES_DIR+'U_b[J].dat')
U_Dcol = lib.read_file_and_define_second_col(RES_DIR+'Ud_col[J].dat')
U_Dres = lib.read_file_and_define_second_col(RES_DIR+'Ud_res[J].dat')
U_EL   = lib.read_file_and_define_second_col(RES_DIR+'U_el[J].dat')
U_MA   = lib.read_file_and_define_second_col(RES_DIR+'U_ma[J].dat')
U_Sd   = lib.read_file_and_define_second_col(RES_DIR+'U_sd[J].dat')
U_Su   = lib.read_file_and_define_second_col(RES_DIR+'U_su[J].dat')
U_Sf   = lib.read_file_and_define_second_col(RES_DIR+'U_sf[J].dat')
U_Sb   = lib.read_file_and_define_second_col(RES_DIR+'U_sb[J].dat')
C_EL   = np.floor(np.log(np.amax(U_B)/np.amax(U_EL))/np.log(10))
C_EL   = 10.**C_EL
C_MA   = np.floor(np.log(np.amax(U_B)/np.amax(U_MA))/np.log(10))
C_MA   = 10.**C_MA
U_stot = []
U_EL_n = []
U_MA_n = []
N_T    = len(T)-1
for i in range(0, N_T+1):
    U_stot.append(U_Su[i]+U_Sd[i]+U_Sf[i]+U_Sb[i])
    U_EL_n.append(C_EL*U_EL[i])
    U_MA_n.append(C_MA*U_MA[i])
balance  = U_E[N_T]-(U_B[N_T]+U_EL[N_T]+U_MA[N_T]+U_Dcol[N_T]+U_Dres[N_T]+U_stot[N_T])
error    = 100.*balance/U_E[N_T]
X_LABEL  = r'$t \, (\mathrm{fs})$'
Y_LABEL  = r'$\mathrm{Energy}\, U\,(\mathrm{J})$'
TTL      = 'Energy conservation error '
TTL     += '= '+str(np.floor(100*error)/100)+r'$\,\%$'
FIG_FILE = SIMU_DIR+'energy_conservation.png'
lib.make_scalars_plot_figure(xplot     = T,
                             yplot1    = U_EL_n,
                             legend1   = r'$U_{E}\times$'+f"{C_EL:.0E}",
                             color1    = 'green',
                             yplot2    = U_MA_n,
                             legend2   = r'$U_{E}\times$'+f"{C_EL:.0E}",
                             color2    = 'cyan',
                             yplot3    = U_E,
                             legend3   = r'$U_\mathrm{inj}$',
                             color3    = 'red',
                             yplot4    = U_B,
                             legend4   = r'$U_b$',
                             color4    = 'black',
                             yplot5    = U_Dcol,
                             legend5   = r'$U_{d,\mathrm{col}}$',
                             color5    = 'magenta',
                             yplot6    = U_Dres,
                             legend6   = r'$U_{d,\mathrm{res}}$',
                             color6    = 'blue',
                             yplot7    = U_stot,
                             legend7   = r'$U_\mathrm{esc}$',
                             color7    = 'red',
                             xplot_min = 0.,
                             xplot_max = T[N_T-1],
                             yplot_min = U_E[0],
                             yplot_max = 1.1*U_E[N_T-1],
                             xlabel    = X_LABEL,
                             ylabel    = Y_LABEL,
                             title     = TTL,
                             filename  = FIG_FILE,
                             logx      = False,
                             logy      = False,
                             grid      = False)
