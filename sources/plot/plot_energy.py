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

T       = lib.read_file_and_define_first_col( RES_DIR+'U_e[J].dat')
U_INJ   = lib.read_file_and_define_second_col(RES_DIR+'U_e[J].dat')
U_B     = lib.read_file_and_define_second_col(RES_DIR+'U_b[J].dat')
U_D_COL = lib.read_file_and_define_second_col(RES_DIR+'Ud_col[J].dat')
U_D_RES = lib.read_file_and_define_second_col(RES_DIR+'Ud_res[J].dat')
U_EL    = lib.read_file_and_define_second_col(RES_DIR+'U_el[J].dat')
U_MA    = lib.read_file_and_define_second_col(RES_DIR+'U_ma[J].dat')
U_S_D   = lib.read_file_and_define_second_col(RES_DIR+'U_sd[J].dat')
U_S_U   = lib.read_file_and_define_second_col(RES_DIR+'U_su[J].dat')
U_S_F   = lib.read_file_and_define_second_col(RES_DIR+'U_sf[J].dat')
U_S_B   = lib.read_file_and_define_second_col(RES_DIR+'U_sb[J].dat')
C_EL    = np.floor(np.log(np.amax(U_B)/np.amax(U_EL))/np.log(10))
C_EL    = 10.**C_EL
C_MA    = np.floor(np.log(np.amax(U_B)/np.amax(U_MA))/np.log(10))
C_MA    = 10.**C_MA
U_stot  = []
U_EL_N  = []
U_MA_N  = []
N_T     = len(T)-1
for i in range(0, N_T+1):
    U_stot.append(U_S_U[i]+U_S_D[i]+U_S_F[i]+U_S_B[i])
    U_EL_N.append(C_EL*U_EL[i])
    U_MA_N.append(C_MA*U_MA[i])
U_BAL    = U_INJ[N_T]-(U_B[N_T]+U_EL[N_T]+U_MA[N_T]+U_D_COL[N_T]+U_D_RES[N_T]+U_stot[N_T])
U_ERR    = 100.*U_BAL/U_INJ[N_T]
T_MIN    = 0.
T_MAX    = T[N_T-1]
U_MIN    = U_INJ[0]
U_MAX    = 1.1*U_INJ[N_T-1]
U_EL_LEG = r'$U_{E}\times$'+f"{C_EL:.0E}"
U_MA_LEG = r'$U_{B}\times$'+f"{C_MA:.0E}"
T_LBL    = r'$t \, (\mathrm{fs})$'
U_LBL    = r'$\mathrm{Energy}\, U\,(\mathrm{J})$'
U_TTL    = 'Energy conservation error '
U_TTL   += '= '+str(np.floor(100*U_ERR)/100)+r'$\,\%$'
U_FIG    = SIMU_DIR+'energy_conservation.png'
lib.make_scalars_plot_figure(xplot     = T,
                             yplot1    = U_EL_N,
                             legend1   = U_EL_LEG,
                             color1    = 'green',
                             yplot2    = U_MA_N,
                             legend2   = U_MA_LEG,
                             color2    = 'cyan',
                             yplot3    = U_INJ,
                             legend3   = r'$U_\mathrm{inj}$',
                             color3    = 'red',
                             yplot4    = U_B,
                             legend4   = r'$U_b$',
                             color4    = 'black',
                             yplot5    = U_D_COL,
                             legend5   = r'$U_{d,\mathrm{col}}$',
                             color5    = 'magenta',
                             yplot6    = U_D_RES,
                             legend6   = r'$U_{d,\mathrm{res}}$',
                             color6    = 'blue',
                             yplot7    = U_stot,
                             legend7   = r'$U_\mathrm{esc}$',
                             color7    = 'red',
                             xplot_min = T_MIN,
                             xplot_max = T_MAX,
                             yplot_min = U_MIN,
                             yplot_max = U_MAX,
                             xlabel    = T_LBL,
                             ylabel    = U_LBL,
                             title     = U_TTL,
                             filename  = U_FIG,
                             logx      = False,
                             logy      = False,
                             grid      = False)
