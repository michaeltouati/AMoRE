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
Read and plot data t (fs) | x (microns) | z (microns) | F from files:
* ne[cm-3].dat
* nb[cm-3].dat
* ni[cm-3].dat
* jb_z[A_cm-2].dat
* jb_x[A_cm-2].dat
* je_z[A_cm-2].dat
* je_x[A_cm-2].dat
* E_x[V_m-1].dat
* E_x[V_m-1].dat
* B_y[Tesla].dat
* We[erg_s-1_cm-3].dat
* Wi[erg_s-1_cm-3].dat
* Te[eV].dat
* Ti[eV].dat
* resis[Ohm.m].dat
* Kappa_e[erg_m-1_K-1_s-1].dat
* K_shell_ioniz_rate_[s-1].dat
* n_Kalpha[cm-3].dat
* n_Kbeta[cm-3].dat
"""
import library as lib

SIMU_NAME=lib.get_results_dir()

lib.create_dir('figures/')

SIMU_DIR = 'figures/'+SIMU_NAME+'/'
lib.create_dir(SIMU_DIR)

MAP_DIR = SIMU_DIR+'2D_maps/'
lib.create_dir(MAP_DIR)

RES_DIR = 'results/'+SIMU_NAME+'/'

print(' ---------------------------------------')
print(' Hydrodynamic Moments 2D Density Plot')
print(' ---------------------------------------')

[N_z,N_x] = lib.find_spatial_simulation_box_dimension(RES_DIR+'B_y[Tesla].dat')
X_LABEL   = r'$z\,(\mu\mathrm{m})$'
Y_LABEL   = r'$x\,(\mu\mathrm{m})$'

plot_ni = lib.get_string_parameter('tabulated_plasma') == '.true.'
plot_ni = plot_ni or  lib.get_string_parameter('Tracer') == '1'
if plot_ni :
    print('  ')
    print(' Target ion density ni')
    NI_RES = RES_DIR+'ni[cm-3].dat'
    NI_TTL = r'$n_i\,(\mathrm{cm}^{-3})$'
    NI_DIR = MAP_DIR+'ni/'
    lib.create_dir(NI_DIR)
    NI_FIG = NI_DIR+'ni_'
    lib.read_and_plot_2d_pcolormesh(res_file = NI_RES,
                                    n_1      = N_x,
                                    n_2      = N_z,
                                    colormap = 'Greens',
                                    title    = NI_TTL,
                                    log      = False,
                                    xlabel   = X_LABEL,
                                    ylabel   = Y_LABEL,
                                    eq_axis  = True,
                                    fig_file = NI_FIG)
print('  ')
print(' Beam density nb')
NB_RES = RES_DIR+'nb[cm-3].dat'
NB_TTL = r'$n_b\,(\mathrm{cm}^{-3})$'
# NB_TTL = r'$\log_{10} ( n_b\,(\mathrm{cm}^{-3}) )$')
NB_DIR = MAP_DIR+'nb/'
lib.create_dir(NB_DIR)
NB_FIG = NB_DIR+'nb_'
lib.read_and_plot_2d_pcolormesh(res_file = NB_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'Blues',
                                title    = NB_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = NB_FIG)
print('  ')
print(' Target electron density ne')
NE_RES = RES_DIR+'ne[cm-3].dat'
NE_TTL = r'$n_e\,(\mathrm{cm}^{-3})$'
NE_DIR = MAP_DIR+'ne/'
lib.create_dir(NE_DIR)
NE_FIG = NE_DIR+'ne_'
lib.read_and_plot_2d_pcolormesh(res_file = NE_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'Reds',
                                title    = NE_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = NE_FIG)
print('  ')
print(' Beam current density jbz')
JBZ_RES = RES_DIR+'jb_z[A_cm-2].dat'
JBZ_TTL = r'$j_{b,z}\,(\mathrm{A/cm}^{2})$'
JBZ_DIR = MAP_DIR+'jbz/'
lib.create_dir(JBZ_DIR)
JBZ_FIG = JBZ_DIR+'jbz_'
lib.read_and_plot_2d_pcolormesh(res_file = JBZ_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = JBZ_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = JBZ_FIG)
print('  ')
print(' Beam current density jbx')
JBX_RES = RES_DIR+'jb_x[A_cm-2].dat'
JBX_TTL = r'$j_{b,x}\,(\mathrm{A/cm}^{2})$'
JBX_DIR = MAP_DIR+'jbx/'
lib.create_dir(JBX_DIR)
JBX_FIG = JBX_DIR+'jbx_'
lib.read_and_plot_2d_pcolormesh(res_file = JBX_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = JBX_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = JBX_FIG)
print('  ')
print(' Return current density jez')
JEZ_RES = RES_DIR+'je_z[A_cm-2].dat'
JEZ_TTL = r'$j_{e,z}\,(\mathrm{A/cm}^{2})$'
JEZ_DIR = MAP_DIR+'jez/'
lib.create_dir(JEZ_DIR)
JEZ_FIG = JEZ_DIR+'jez_'
lib.read_and_plot_2d_pcolormesh(res_file = JEZ_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = JEZ_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = JEZ_FIG)
print('  ')
print(' Return current density jex')
JEX_RES = RES_DIR+'je_x[A_cm-2].dat'
JEX_TTL = r'$j_{e,x}\,(\mathrm{A/cm}^{2})$'
JEX_DIR = MAP_DIR+'jex/'
lib.create_dir(JEX_DIR)
JEX_FIG = JEX_DIR+'jex_'
lib.read_and_plot_2d_pcolormesh(res_file = JEX_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = JEX_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = JEX_FIG)
print('  ')
print(' Electric field Ez')
EZ_RES = RES_DIR+'E_z[V_m-1].dat'
EZ_TTL = r'$E_{z}\,(\mathrm{V/m})$'
EZ_DIR = MAP_DIR+'Ez/'
lib.create_dir(EZ_DIR)
EZ_FIG=EZ_DIR+'Ez_'
lib.read_and_plot_2d_pcolormesh(res_file = EZ_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = EZ_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = EZ_FIG)
print('  ')
print(' Electric field Ex')
EX_RES = RES_DIR+'E_x[V_m-1].dat'
EX_TTL = r'$E_{x}\,(\mathrm{V/m})$'
EX_DIR = MAP_DIR+'Ex/'
lib.create_dir(EX_DIR)
EX_FIG = EX_DIR+'Ex_'
lib.read_and_plot_2d_pcolormesh(res_file = EX_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = EX_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = EX_FIG)
print('  ')
print(' Magnetic field By')
BY_RES = RES_DIR+'B_y[Tesla].dat'
BY_TTL = r'$B_y\,(\mathrm{T})$'
BY_DIR = MAP_DIR+'By/'
lib.create_dir(BY_DIR)
BY_FIG = BY_DIR+'By_'
lib.read_and_plot_2d_pcolormesh(res_file = BY_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'seismic',
                                title    = BY_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = BY_FIG)
print('  ')
print(' Density of power deposited on background electrons We')
WE_RES = RES_DIR+'We[erg_s-1_cm-3].dat'
WE_TTL = r'$W_e\,(\mathrm{erg/s.cm}^{3})$'
WE_DIR = MAP_DIR+'We/'
lib.create_dir(WE_DIR)
WE_FIG = WE_DIR+'We_'
lib.read_and_plot_2d_pcolormesh(res_file = WE_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'OrRd',
                                title    = WE_TTL,
                                log      = False,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = WE_FIG)
if lib.get_string_parameter('bi_temp') == '.true.':
    print('  ')
    print(' Density of power deposited on background ions Wi')
    WI_RES = RES_DIR+'Wi[erg_s-1_cm-3].dat'
    WI_TTL = r'$W_i\,(\mathrm{erg/s.cm}^{3})$'
    WI_DIR = MAP_DIR+'Wi/'
    lib.create_dir(WI_DIR)
    WI_FIG = WI_DIR+'Wi_'
    lib.read_and_plot_2d_pcolormesh(res_file = WI_RES,
                                    n_1      = N_x,
                                    n_2      = N_z,
                                    colormap = 'BuGn',
                                    title    = WI_TTL,
                                    log      = False,
                                    xlabel   = X_LABEL,
                                    ylabel   = Y_LABEL,
                                    eq_axis  = True,
                                    fig_file = WI_FIG)
print('  ')
print(' Background electron temperature Te')
TE_RES = RES_DIR+'Te[eV].dat'
TE_TTL = r'$\log_{10}{\left(T_e\,(\mathrm{eV})\right)}$'
TE_DIR = MAP_DIR+'Te/'
lib.create_dir(TE_DIR)
TE_FIG = TE_DIR+'Te_'
lib.read_and_plot_2d_pcolormesh(res_file = TE_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'jet',
                                title    = TE_TTL,
                                log      = True,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = TE_FIG)
if lib.get_string_parameter('bi_temp') == '.true.':
    print('  ')
    print(' Background ion temperature Ti')
    TI_RES = RES_DIR+'Ti[eV].dat'
    TI_TTL = r'$\log_{10}{\left(T_i\,(\mathrm{eV})\right)}$'
    TI_DIR = MAP_DIR+'Ti/'
    lib.create_dir(TI_DIR)
    TI_FIG = TI_DIR+'Ti_'
    lib.read_and_plot_2d_pcolormesh(res_file = TI_RES,
                                    n_1      = N_x,
                                    n_2      = N_z,
                                    colormap = 'jet',
                                    title    = TI_TTL,
                                    log      = True,
                                    xlabel   = X_LABEL,
                                    ylabel   = Y_LABEL,
                                    eq_axis  = True,
                                    fig_file = TI_FIG)
print('  ')
print(' Electrical resistivity eta')
ETA_RES = RES_DIR+'resis[Ohm.m].dat'
ETA_TTL = r'$\log_{10}{\left (\eta\,(\Omega.\mathrm{m}) \right )}$'
ETA_DIR = MAP_DIR+'eta/'
lib.create_dir(ETA_DIR)
ETA_FIG = ETA_DIR+'eta_'
lib.read_and_plot_2d_pcolormesh(res_file = ETA_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'terrain',
                                title    = ETA_TTL,
                                log      = True,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = ETA_FIG)
print('  ')
print(' Thermal conductivity kappa')
KAPPA_RES  = RES_DIR+'Kappa_e[erg_m-1_K-1_s-1].dat'
KAPPA_TTL  = r'$\log_{10}\Big (\kappa\,$'
KAPPA_TTL += r'$(\mathrm{erg/m.K.s}) \Big )$'
KAPPA_DIR  = MAP_DIR+'kappa/'
lib.create_dir(KAPPA_DIR)
KAPPA_FIG  = KAPPA_DIR+'kappa_'
lib.read_and_plot_2d_pcolormesh(res_file = KAPPA_RES,
                                n_1      = N_x,
                                n_2      = N_z,
                                colormap = 'terrain',
                                title    = KAPPA_TTL,
                                log      = True,
                                xlabel   = X_LABEL,
                                ylabel   = Y_LABEL,
                                eq_axis  = True,
                                fig_file = KAPPA_FIG)
if lib.get_string_parameter('Kalpha') == '.true.':
    print('  ')
    print(' Ionization rate of K-shell electrons nuK')
    NUK_RES = RES_DIR+'K_shell_ioniz_rate_[s-1].dat'
    NUK_TTL = r'$\nu_K\,(\mathrm{s}^{-1})$'
    NUK_DIR = MAP_DIR+'nuK/'
    lib.create_dir(NUK_DIR)
    NUK_FIG = NUK_DIR+'nuK_'
    lib.read_and_plot_2d_pcolormesh(res_file = NUK_RES,
                                    n_1      = N_x,
                                    n_2      = N_z,
                                    colormap = 'OrRd',
                                    title    = NUK_TTL,
                                    log      = False,
                                    xlabel   = X_LABEL,
                                    ylabel   = Y_LABEL,
                                    eq_axis  = True,
                                    fig_file = NUK_FIG)
    print('  ')
    print(' Time integrated density of emitted Kalpha photons nKalpha')
    NKA_RES = RES_DIR+'n_Kalpha[cm-3].dat'
    # NKA_TTL = r'$n_{K\alpha}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1})$'
    NKA_TTL = r'$\log_{10} ( n_{K\alpha}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1}))$'
    NKA_DIR = MAP_DIR+'nKa/'
    lib.create_dir(NKA_DIR)
    NKA_FIG = NKA_DIR+'nKa_'
    lib.read_and_plot_2d_pcolormesh(res_file = NKA_RES,
                                    n_1      = N_x,
                                    n_2      = N_z,
                                    colormap = 'hot',
                                    title    = NKA_TTL,
                                    log      = True,
                                    xlabel   = X_LABEL,
                                    ylabel   = Y_LABEL,
                                    eq_axis  = True,
                                    fig_file = NKA_FIG)
    print('  ')
    print(' Time integrated density of emitted Kbeta photons nKb')
    NKB_RES = RES_DIR+'n_Kbeta[cm-3].dat'
    # NKB_TTL = r'$n_{K\beta}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1})$'
    NKB_TTL = r'$\log_{10} ( n_{K\beta}\,(\mathrm{cm}^{-3}.\mathrm{sr}^{-1}) )$'
    NKB_DIR = MAP_DIR+'nKb/'
    lib.create_dir(NKB_DIR)
    NKB_FIG = NKB_DIR+'nKb_'
    lib.read_and_plot_2d_pcolormesh(res_file = NKB_RES,
                                    n_1      = N_x,
                                    n_2      = N_z,
                                    colormap = 'hot',
                                    title    = NKB_TTL,
                                    log      = True,
                                    xlabel   = X_LABEL,
                                    ylabel   = Y_LABEL,
                                    eq_axis  = True,
                                    fig_file = NKB_FIG)
