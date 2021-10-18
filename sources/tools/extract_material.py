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
import rapc

dir= os.path.dirname("figures/")
if not os.path.exists(dir):
    os.mkdir(dir)

dir= os.path.dirname("figures/material_properties/")
if not os.path.exists(dir):
    os.mkdir(dir)
    
rapc.read_and_plot_curve('diag/zeffvsTe[eV].dat',
							r'$\mathbf{T_e}\,(\mathrm{eV})$',
							r'$\mathbf{Z^*}\,()$',
							'Ionization state',
							'figures/material_properties/Zf.png',
							1,0)
rapc.read_and_plot_curve('diag/stopping_power[keV_Microns-1]vsEps[keV].dat',
							r'$\mathbf{\varepsilon}\,(\mathrm{keV})$',
							r'$\mathbf{S}\,(\mathrm{keV}/\mu\mathrm{m})$',
							'Beam electron stopping power',
							'figures/material_properties/S.png',
							1,1)
rapc.read_and_plot_curve('diag/ang_coll_rate[s-1]vsEps[keV].dat',
							r'$\mathbf{\varepsilon}\,(\mathrm{keV})$',
							r'$\mathbf{\nu}\,(\mathrm{s}^{-1})$',
							'Beam electron isotropization rate',
							'figures/material_properties/nu.png',
							1,1)
rapc.read_and_plot_two_log_curve('diag/resistivity[SI]vsTe[eV]_Te_eq_Ti.dat',
							'diag/resistivity[SI]vsTe[eV]_Ti_eq_Tamb.dat',
							r'$\mathbf{T_i}=\mathbf{T_e}$',
							r'$\mathbf{T_i}=300\,\mathrm{K}$',
							'upper right',
							r'$\mathbf{T_e}$',
							r'$\mathbf{\eta}\,(\mathrm{\Omega.m})$',
							'Electrical resistivity',
							'figures/material_properties/eta.png')
rapc.read_and_plot_two_log_curve('diag/conductivity[SI]vsTe[eV]_Te_eq_Ti.dat',
							'diag/conductivity[SI]vsTe[eV]_Ti_eq_Tamb.dat',
							r'$\mathbf{T_i}=\mathbf{T_e}$',
							r'$\mathbf{T_i}=300\,\mathrm{K}$',
							'upper left',
							r'$\mathbf{T_e}\,(\mathrm{eV})$',
							r'$\mathbf{\kappa}\,(\mathrm{J/m/K/s})$',
							'Thermal conductivity',
							'figures/material_properties/kappa.png')
rapc.read_and_plot_two_log_curve('diag/G[SI]vsTe[eV]_Te_eq_Ti.dat',
							'diag/G[SI]vsTe[eV]_Ti_eq_Tamb.dat',
							r'$\mathbf{T_i}=\mathbf{T_e}$',
							r'$\mathbf{T_i}=300\,\mathrm{K}$',
							'lower center',
							r'$\mathbf{T_e}\,(\mathrm{eV})$',
							r'$\mathbf{\Omega_{ei}}\,(\mathrm{J}/\mathrm{m}^3/\mathrm{s/K})$',
							'Electron-ion coupling factor',
							'figures/material_properties/omega_ei.png')
rapc.read_and_plot_two_log_curve('diag/electron_capacity[SI]vsTe[eV].dat',
							'diag/ion_capacity[SI]vsTe[eV].dat',
							r'$\mathbf{C_{V,e}}$',
							r'$\mathbf{C_{V,i}}$',
							'upper left',
							r'$\mathbf{T=T_e}\,\mathrm{or}\,\mathbf{T_i}\,(\mathrm{eV})$',
							r'$\mathbf{C_V}\,(\mathrm{J}/\mathrm{m}^3.\mathrm{K})$',
							'Thermal capacities',
							'figures/material_properties/Cv.png')