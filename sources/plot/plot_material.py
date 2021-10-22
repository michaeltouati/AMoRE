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
import library as lib

simu_name=lib.get_results_dir()

print(' -------------------------------')
print(' Material Properties Scalar Plot')
print(' -------------------------------')
print('  ')

lib.create_dir('figures/')
lib.create_dir('figures/'+simu_name+'/')

subdir = 'figures/'+simu_name+'/material_properties/'
lib.create_dir(subdir)

results_dir = 'results/'+simu_name+'/'

print(' Ionization state Z*')
lib.read_and_plot_curve(results_dir+'zeffvsTe[eV].dat',
						r'$\mathbf{T_e}\,(\mathrm{eV})$',
						r'$\mathbf{Z^*}\,()$',
						'Ionization state',
						subdir+'Zf.png',
						1,0)
print('  ')
print(' Stopping power S')
lib.read_and_plot_curve(results_dir+'stopping_power[keV_Microns-1]vsEps[keV].dat',
						r'$\mathbf{\varepsilon}\,(\mathrm{keV})$',
						r'$\mathbf{S}\,(\mathrm{keV}/\mu\mathrm{m})$',
						'Beam electron stopping power',
						subdir+'S.png',
						1,1)
print('  ')
print(' Angular isotropization rate nu')
lib.read_and_plot_curve(results_dir+'ang_coll_rate[s-1]vsEps[keV].dat',
						r'$\mathbf{\varepsilon}\,(\mathrm{keV})$',
						r'$\mathbf{\nu}\,(\mathrm{s}^{-1})$',
						'Beam electron isotropization rate',
						subdir+'nu.png',
						1,1)
print('  ')
print(' Electrical resistivity eta')
lib.read_and_plot_two_log_curve(results_dir+'resistivity[SI]vsTe[eV]_Te_eq_Ti.dat',
								results_dir+'resistivity[SI]vsTe[eV]_Ti_eq_Tamb.dat',
								r'$\mathbf{T_i}=\mathbf{T_e}$',
								r'$\mathbf{T_i}=293.15\,\mathrm{K}$',
								'upper right',
								r'$\mathbf{T_e}$',
								r'$\mathbf{\eta}\,(\mathrm{\Omega.m})$',
								'Electrical resistivity',
								subdir+'eta.png')
print('  ')
print(' Thermal conductivity kappa')
lib.read_and_plot_two_log_curve(results_dir+'conductivity[SI]vsTe[eV]_Te_eq_Ti.dat',
								results_dir+'conductivity[SI]vsTe[eV]_Ti_eq_Tamb.dat',
								r'$\mathbf{T_i}=\mathbf{T_e}$',
								r'$\mathbf{T_i}=293.15\,\mathrm{K}$',
								'upper left',
								r'$\mathbf{T_e}\,(\mathrm{eV})$',
								r'$\mathbf{\kappa}\,(\mathrm{J/m/K/s})$',
								'Thermal conductivity',
								subdir+'kappa.png')
print('  ')
print(' Electron-ion/lattice coupling factor G')
lib.read_and_plot_two_log_curve(results_dir+'G[SI]vsTe[eV]_Te_eq_Ti.dat',
								results_dir+'G[SI]vsTe[eV]_Ti_eq_Tamb.dat',
								r'$\mathbf{T_i}=\mathbf{T_e}$',
								r'$\mathbf{T_i}=293.15\,\mathrm{K}$',
								'lower center',
								r'$\mathbf{T_e}\,(\mathrm{eV})$',
								r'$\mathbf{\Omega_{ei}}\,(\mathrm{J}/\mathrm{m}^3/\mathrm{s/K})$',
								'Electron-ion coupling factor',
								subdir+'omega_ei.png')
print('  ')
print(' Target electrons and ions thermal capacities Cve and Cvi')
lib.read_and_plot_two_log_curve(results_dir+'electron_capacity[SI]vsTe[eV].dat',
								results_dir+'ion_capacity[SI]vsTe[eV].dat',
								r'$\mathbf{C_{V,e}}$',
								r'$\mathbf{C_{V,i}}$',
								'upper left',
								r'$\mathbf{T=T_e}\,\mathrm{or}\,\mathbf{T_i}\,(\mathrm{eV})$',
								r'$\mathbf{C_V}\,(\mathrm{J}/\mathrm{m}^3.\mathrm{K})$',
								'Thermal capacities',
								subdir+'Cv.png')
print('  ')