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
Read and plot data :
* E (kEV) | dE_dS (keV/microns) / nu (/s)
  * stopping_power[keV_Microns-1]vsEps[keV].dat
  * ang_coll_rate[s-1]vsEps[keV].dat
* dT (eV) | Plasma transport coefficient (SI) from files :
  * zeffvsTe[eV].dat
  * resistivity[SI]vsTe[eV]_Te_eq_Ti.dat
  * resistivity[SI]vsTe[eV]_Ti_eq_Tamb.dat
  * conductivity[SI]vsTe[eV]_Te_eq_Ti.dat
  * conductivity[SI]vsTe[eV]_Ti_eq_Tamb.dat
  * G[SI]vsTe[eV]_Te_eq_Ti.dat
  * G[SI]vsTe[eV]_Ti_eq_Tamb.dat
  * electron_capacity[SI]vsTe[eV].dat
  * ion_capacity[SI]vsTe[eV].dat
"""
import library as lib
import numpy as np

SIMU_NAME=lib.get_results_dir()

print(' -------------------------------')
print(' Material Properties Scalar Plot')
print(' -------------------------------')
print('  ')

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

MAT_DIR = 'figures/'+SIMU_NAME+'/material_properties/'
lib.create_dir(MAT_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

print(' Stopping power dE/ds')
print('  ')
RES_FILE      = RES_DIR+'stopping_power[keV_Microns-1]vsEps[keV].dat'
[EPS,DEPS_DS] = lib.read_file_and_define_two_first_cols(RES_FILE)
EPS_MIN       = lib.get_log_axis_min_value(EPS)
EPS_MAX       = lib.get_log_axis_max_value(EPS)
DEPS_DS_MIN   = lib.get_log_axis_min_value(DEPS_DS)
DEPS_DS_MAX   = lib.get_log_axis_max_value(DEPS_DS)
FIG_FILE      = MAT_DIR+'dE_ds.png'
lib.make_one_scalar_plot_figure(xplot=EPS, yplot=DEPS_DS,
	                              color='red',
	                              xplot_min=EPS_MIN,
	                              xplot_max=EPS_MAX,
	                              yplot_min=DEPS_DS_MIN,
	                              yplot_max=DEPS_DS_MAX,
	                              xlabel=r'$\varepsilon\,(\mathrm{keV})$',
	                              ylabel=r'$d \varepsilon / d s \,(\mathrm{keV}/\mu\mathrm{m})$',
	                              title='Beam electrons stopping power',
	                              filename=FIG_FILE,
	                              logx=True, logy=True, grid=True)

print(' Angular isotropization rate nu')
print('  ')
RES_FILE = RES_DIR+'ang_coll_rate[s-1]vsEps[keV].dat'
[EPS,NU] = lib.read_file_and_define_two_first_cols(RES_FILE)
EPS_MIN  = lib.get_log_axis_min_value(EPS)
EPS_MAX  = lib.get_log_axis_max_value(EPS)
NU_MIN   = lib.get_log_axis_min_value(NU)
NU_MAX   = lib.get_log_axis_max_value(NU)
EPS_LBL  = r'$\varepsilon\,(\mathrm{keV})$'
NU_LBL   = r'$\nu\,(\mathrm{s}^{-1})$'
TTL      = 'Beam electrons isotropization rate'
FIG_FILE = MAT_DIR+'nu.png'
lib.make_one_scalar_plot_figure(xplot=EPS, yplot=NU,
	                              color='red',
	                              xplot_min=EPS_MIN,
	                              xplot_max=EPS_MAX,
	                              yplot_min=NU_MIN,
	                              yplot_max=NU_MAX,
	                              xlabel=EPS_LBL,
	                              ylabel=NU_LBL,
	                              title=TTL,
	                              filename=FIG_FILE,
	                              logx=True, logy=True, grid=True)

print(' Ionization state Z*')
print('  ')
RES_FILE = RES_DIR+'zeffvsTe[eV].dat'
[TE,ZF]  = lib.read_file_and_define_two_first_cols(RES_FILE)
TE_MIN   = lib.get_log_axis_min_value(TE)
TE_MAX   = lib.get_log_axis_max_value(TE)
ZF_MIN   = 0.
ZF_MAX   = np.ceil(np.max(ZF))
TE_LBL   = r'$T_e\,(\mathrm{eV})$'
ZF_LBL   = r'$Z^*\,(\,)$'
TTL      = 'Ionization state'
FIG_FILE = MAT_DIR+'Zf.png'
lib.make_one_scalar_plot_figure(xplot=TE, yplot=ZF,
	                              color='red',
	                              xplot_min=TE_MIN,
	                              xplot_max=TE_MAX,
	                              yplot_min=ZF_MIN,
	                              yplot_max=ZF_MAX,
	                              xlabel=TE_LBL,
	                              ylabel=ZF_LBL,
	                              title=TTL,
	                              filename=FIG_FILE,
	                              logx=True, logy=False, grid=True)

print(' Electrical resistivity eta')
print('  ')
RES_FILE1 = RES_DIR+'resistivity[SI]vsTe[eV]_Te_eq_Ti.dat'
RES_FILE2 = RES_DIR+'resistivity[SI]vsTe[eV]_Ti_eq_Tamb.dat'
[TE,ETA1] = lib.read_file_and_define_two_first_cols(RES_FILE1)
[TE,ETA2] = lib.read_file_and_define_two_first_cols(RES_FILE2)
TE_MIN    = lib.get_log_axis_min_value(TE)
TE_MAX    = lib.get_log_axis_max_value(TE)
ETA1_MIN  = lib.get_log_axis_min_value(ETA1)
ETA1_MAX  = lib.get_log_axis_max_value(ETA1)
ETA2_MIN  = lib.get_log_axis_min_value(ETA2)
ETA2_MAX  = lib.get_log_axis_max_value(ETA2)
ETA_MIN   = min(ETA1_MIN,ETA2_MIN)
ETA_MAX   = max(ETA1_MAX,ETA2_MAX)
TE_LBL    = r'$T_e\,(\mathrm{eV})$'
ETA_LBL   = r'$\eta\,(\mathrm{\Omega.m})$'
TTL       = 'Electrical resistivity'
FIG_FILE  = MAT_DIR+'eta.png'
lib.make_scalars_plot_figure(xplot=TE,
	                           yplot1=ETA1,
	                           legend1 =r'$T_i=T_e$',
	                           color1='red',
	                           yplot2=ETA2,
	                           legend2 =r'$T_i=293.15\,\mathrm{K}$',
	                           color2='blue',
	                           xplot_min=TE_MIN,
	                           xplot_max=TE_MAX,
	                           yplot_min=ETA_MIN,
	                           yplot_max=ETA_MAX,
	                           xlabel=TE_LBL,
	                           ylabel=ETA_LBL,
	                           title=TTL,
	                           filename=FIG_FILE,
	                           logx=True, logy=True, grid=True)

print(' Thermal conductivity kappa')
print('  ')
RES_FILE1   = RES_DIR+'conductivity[SI]vsTe[eV]_Te_eq_Ti.dat'
RES_FILE2   = RES_DIR+'conductivity[SI]vsTe[eV]_Ti_eq_Tamb.dat'
[TE,KAPPA1] = lib.read_file_and_define_two_first_cols(RES_FILE1)
[TE,KAPPA2] = lib.read_file_and_define_two_first_cols(RES_FILE2)
TE_MIN      = lib.get_log_axis_min_value(TE)
TE_MAX      = lib.get_log_axis_max_value(TE)
KAPPA1_MIN  = lib.get_log_axis_min_value(KAPPA1)
KAPPA1_MAX  = lib.get_log_axis_max_value(KAPPA1)
KAPPA2_MIN  = lib.get_log_axis_min_value(KAPPA2)
KAPPA2_MAX  = lib.get_log_axis_max_value(KAPPA2)
KAPPA_MIN   = min(KAPPA1_MIN,KAPPA2_MIN)
KAPPA_MAX   = max(KAPPA1_MAX,KAPPA2_MAX)
TE_LBL      = r'$T_e\,(\mathrm{eV})$'
KAPPA_LBL   = r'$\kappa_e\,(\mathrm{J/m/K/s})$'
TTL         = 'Electron thermal conductivity'
FIG_FILE    = MAT_DIR+'kappa.png'
lib.make_scalars_plot_figure(xplot=TE,
	                           yplot1=KAPPA1,
	                           legend1 =r'$T_i=T_e$',
	                           color1='red',
	                           yplot2=KAPPA2,
	                           legend2 =r'$T_i=293.15\,\mathrm{K}$',
	                           color2='blue',
	                           xplot_min=TE_MIN,
	                           xplot_max=TE_MAX,
	                           yplot_min=KAPPA_MIN,
	                           yplot_max=KAPPA_MAX,
	                           xlabel=TE_LBL,
	                           ylabel=KAPPA_LBL,
	                           title=TTL,
	                           filename=FIG_FILE,
	                           logx=True, logy=True, grid=True)

print(' Electron-ion/lattice coupling factor G')
print('  ')
RES_FILE1 = RES_DIR+'G[SI]vsTe[eV]_Te_eq_Ti.dat'
RES_FILE2 = RES_DIR+'G[SI]vsTe[eV]_Ti_eq_Tamb.dat'
[TE,G1]   = lib.read_file_and_define_two_first_cols(RES_FILE1)
[TE,G2]   = lib.read_file_and_define_two_first_cols(RES_FILE2)
TE_MIN    = lib.get_log_axis_min_value(TE)
TE_MAX    = lib.get_log_axis_max_value(TE)
G1_MIN    = lib.get_log_axis_min_value(G1)
G1_MAX    = lib.get_log_axis_max_value(G1)
G2_MIN    = lib.get_log_axis_min_value(G2)
G2_MAX    = lib.get_log_axis_max_value(G2)
G_MIN     = min(G1_MIN,G2_MIN)
G_MAX     = max(G1_MAX,G2_MAX)
TE_LBL    = r'$T_e\,(\mathrm{eV})$'
G_LBL     = r'$\Omega_{ei}\,(\mathrm{J}/\mathrm{m}^3/\mathrm{s/K})$'
TTL       = 'Electron-ion/lattice coupling factor'
FIG_FILE  = MAT_DIR+'Omega_ei.png'
lib.make_scalars_plot_figure(xplot=TE,
	                           yplot1=G1,
	                           legend1 =r'$T_i=T_e$',
	                           color1='red',
	                           yplot2=G2,
	                           legend2 =r'$T_i=293.15\,\mathrm{K}$',
	                           color2='blue',
	                           xplot_min=TE_MIN,
	                           xplot_max=TE_MAX,
	                           yplot_min=G_MIN,
	                           yplot_max=G_MAX,
	                           xlabel=TE_LBL,
	                           ylabel=G_LBL,
	                           title=TTL,
	                           filename=FIG_FILE,
	                           logx=True, logy=True, grid=True)

print(' Target electrons and ions thermal capacities Cve and Cvi')
print('  ')
RES_FILE_E = RES_DIR+'electron_capacity[SI]vsTe[eV].dat'
RES_FILE_I = RES_DIR+'ion_capacity[SI]vsTi[eV].dat'
[TE,CV_E]  = lib.read_file_and_define_two_first_cols(RES_FILE_E)
[TE,CV_I]  = lib.read_file_and_define_two_first_cols(RES_FILE_I)
TE_MIN     = lib.get_log_axis_min_value(TE)
TE_MAX     = lib.get_log_axis_max_value(TE)
CV_E_MIN   = lib.get_log_axis_min_value(CV_E)
CV_E_MAX   = lib.get_log_axis_max_value(CV_E)
CV_I_MIN   = lib.get_log_axis_min_value(CV_I)
CV_I_MAX   = lib.get_log_axis_max_value(CV_I)
CV_MIN     = min(CV_E_MIN,CV_I_MIN)
CV_MAX     = max(CV_E_MAX,CV_I_MAX)
TE_LBL     = r'$T_{e,i}\,(\mathrm{eV})$'
CV_LBL     = r'$C_{V_{e,i}}\,(\mathrm{J}/\mathrm{m}^3.\mathrm{K})$'
TTL        = 'Thermal capacities'
FIG_FILE   = MAT_DIR+'Cv.png'
lib.make_scalars_plot_figure(xplot=TE,
	                           yplot1=CV_E,
	                           legend1 =r'$C_{V_e}$',
	                           color1='red',
	                           yplot2=CV_I,
	                           legend2 =r'$C_{V_i}$',
	                           color2='blue',
	                           xplot_min=TE_MIN,
	                           xplot_max=TE_MAX,
	                           yplot_min=CV_MIN,
	                           yplot_max=CV_MAX,
	                           xlabel=TE_LBL,
	                           ylabel=CV_LBL,
	                           title=TTL,
	                           filename=FIG_FILE,
	                           logx=True, logy=True, grid=True)
