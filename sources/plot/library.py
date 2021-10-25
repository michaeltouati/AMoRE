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
from matplotlib import cm
import os

def get_results_dir():
    file = open('input-deck','r')
    for line in file:
        line      = line.strip()
        array     = line.split()
        if (array[0] == '#simu') :
            name = array[1]
            break
    file.close()
    to_print=' '+name+' AMoRE SIMULATION PLOTS :'
    N_to_print = len(to_print)
    line_to_print= ' '
    for char in range(0,len(to_print)-1):
        line_to_print+='='
    print(line_to_print)
    print(to_print)
    print(line_to_print)
    return name

def create_dir(name):
    dir_name= os.path.dirname(name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def find_spatial_simulation_box_dimension(file_name):
	t0  = []
	x0  = []
	nb_file0 = open(file_name, 'r')
	line = nb_file0.readline()
	line = line.strip()
	array = line.split()
	t0.append(float(array[0]))
	x0.append(float(array[1]))
	counter = 0
	for line in nb_file0:
		line      = line.strip()
		array     = line.split()
		t0.append(float(array[0]))
		x0.append(float(array[1]))
		counter = counter + 1
		if x0[counter]!=x0[counter-1]:
			Nz = counter
			break
	for line in nb_file0:
		line      = line.strip()
		array     = line.split()
		t0.append(float(array[0]))
		counter = counter + 1
		if t0[counter]!=t0[counter-1]:
			NxNz = counter
			break
	nb_file0.close()	
	Nx = NxNz / Nz
	return [Nz,Nx]

def find_energy_bins(file_name):
	x0  = []
	e0  = []
	psi0_file0 = open(file_name, 'r')
	line = psi0_file0.readline()
	line = line.strip()
	array = line.split()
	x0.append(float(array[1]))
	e0.append(float(array[2]))
	counter = 0
	for line in psi0_file0:
		line      = line.strip()
		array     = line.split()
		x0.append(float(array[1]))
		e0.append(float(array[2]))
		counter = counter + 1
		if x0[counter]!=x0[counter-1]:
			Neps = counter
			break
	psi0_file0.close()
	return [Neps,e0]

def read_and_plot_curve(filename,xlabel,ylabel,title,name,logx,logy,font_size,font):
	[x,y] = read_file_and_define_two_first_cols(filename) 
	fig=plt.figure()
	plt.rc('text', usetex=True)
	if logx==1:
		if logy ==0:
			plt.semilogx(x, y,linewidth=2)	
		else:
			plt.loglog(x, y,linewidth=2)
	plt.title(title, fontdict=font)
	plt.xticks(fontsize=font_size)
	plt.xlabel(xlabel, fontdict=font)
	plt.ylabel(ylabel, fontdict=font)
	plt.yticks(fontsize=font_size)
	plt.grid(which='both', axis='both')
	fig.savefig(name,bbox_inches='tight')
	plt.close(fig)
	return;

def read_and_plot_two_log_curve(filename1,filename2,label1,label2,loc,xlabel,ylabel,title,name,font_size,font):
	[x1,y1] = read_file_and_define_two_first_cols(filename1)
	[x2,y2] = read_file_and_define_two_first_cols(filename2)
	fig=plt.figure()
	plt.rc('text', usetex=True)
	plt.loglog(x1, y1,'red',linewidth=2,label=label1)
	plt.loglog(x2, y2,'blue',linewidth=2,label=label2)
	plt.title(title, fontdict=font)
	plt.xticks(fontsize=font_size)
	plt.xlabel(xlabel, fontdict=font)
	plt.ylabel(ylabel, fontdict=font)
	plt.yticks(fontsize=font_size)
	leg = plt.legend(loc=loc,fontsize=font_size, fancybox=True)
	leg.get_frame().set_alpha(0.5)
	plt.grid(which='both', axis='both')
	fig.savefig(name,bbox_inches='tight')
	plt.close(fig)
	return;

def read_file_and_define_first_col(filename,time):
	file = open(filename, 'r')
	for line in file:
		line      = line.strip()
		array     = line.split()
		time.append(float(array[0]))
	file.close()   
	return;

def read_file_and_define_second_col(filename,vector):
	file = open(filename, 'r')
	for line in file:
		line      = line.strip()
		array     = line.split()
		vector.append(float(array[1]))
	file.close()   
	return;

def read_file_and_define_two_first_cols(filename):
	x0 = []
	y0 = []
	file = open(filename, 'r')
	for line in file:
		line      = line.strip()
		array     = line.split()
		x0.append(float(array[0]))
		y0.append(float(array[1]))
	file.close()
	return [x0,y0]

def read_and_plot_2D_pcolormesh(filename,N1,N2,cmap,title,name,log):
	N1 = int(N1)
	N2 = int(N2)
	N3 = N1*N2
	x = []
	z = []
	p = []
	Z = np.zeros((N1,N2))
	X = np.zeros((N1,N2))
	P = np.zeros((N1,N2))
	file = open(filename, 'r')
	counter = 0
	for line in file:
		line      = line.strip()
		array     = line.split()
		x.append(float(array[1]))
		z.append(float(array[2]))
		if log==1:
			p.append(np.log(float(array[3]))/np.log(10))
		else:
			p.append(float(array[3]))
		counter = counter + 1
		if counter % N3 == 0:
			time = int(100.*float(array[0]))/100.
			print(' * t = '+str(time)+' fs')
			N = int(counter / N3)
			val=max(abs(np.amax(p[(N-1)*N3:N*N3])),abs(np.amin(p[(N-1)*N3:N*N3])))
			if (val != 0.) and (log!=1):
				Power=math.floor(np.log(val)/np.log(10.))
			else:
				Power = 0.
			p[(N-1)*N3:N*N3] = [pp / (10.**Power) for pp in p[(N-1)*N3:N*N3]]
			for i in range(0,N1):
				for k in range(0,N2):
					Z[i][k]=z[(N-1)*N3+i*N2+k]
					X[i][k]=x[(N-1)*N3+i*N2+k]
					P[i][k]=p[(N-1)*N3+i*N2+k]
			cmap = plt.get_cmap(cmap)
			if filename in ('diag/jb_z[A_cm-2].dat','diag/je_z[A_cm-2].dat','diag/E_z[V_m-1].dat','diag/B_y[Tesla].dat'):
				Maxval = max(np.amax(p[(N-1)*N3:N*N3]),abs(np.amin(p[(N-1)*N3:N*N3])))
				Minval = -Maxval
			else:
				Maxval = np.amax(p[(N-1)*N3:N*N3])
				Minval = np.amin(p[(N-1)*N3:N*N3])
			norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
			fig=plt.figure()
			plt.rc('text', usetex=True)
			plt.pcolormesh(Z,X,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
			cbar=plt.colorbar()
			cbar.ax.tick_params(labelsize=16)
			plt.gca().set_aspect('equal')
			#font = {'family': 'non-serif',
			font = {'style':  'normal',
        			'color':  'black',
       		 		'weight': 'normal',
       		 		'size': 16,
       		 	}
			if log!=1:
				plt.title(str('%.0e' % (10.**(-Power)))+r'$\times$'+title+' at '+str(time)+' fs', fontdict=font)
			else:
				plt.title(title+' at '+str(time)+' fs', fontdict=font)
			plt.xticks(fontsize=16)
			plt.xlabel(r'$z\,(\mu\mathrm{m})$', fontdict=font)
			plt.xlim([np.amin(z[(N-1)*N3:N*N3]),np.amax(z[(N-1)*N3:N*N3])])
			plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=font)
			plt.yticks(fontsize=16)
			plt.ylim([np.amin(x[(N-1)*N3:N*N3]),np.amax(x[(N-1)*N3:N*N3])])
			fig.savefig(name+str(N)+'.png',bbox_inches='tight')
			plt.close(fig)
	file.close()
	return;
	
def read_and_plot_2D_pcolormesh_abs(filenamez,filenamex,N1,N2,cmap,title,name):
	N1 = int(N1)
	N2 = int(N2)
	N3 = N1*N2
	x = []
	z = []
	pz = []
	px = []
	Z = np.zeros((N1,N2))
	X = np.zeros((N1,N2))
	p  = np.zeros((N1,N2))
	P  = np.zeros((N1,N2))
	filez = open(filenamez, 'r')
	filex = open(filenamex, 'r')
	counter = 0
	for linez in filez:
		linez      = linez.strip()
		arrayz     = linez.split()
		x.append(float(arrayz[1]))
		z.append(float(arrayz[2]))
		pz.append(float(arrayz[3]))
		counter = counter + 1
		if counter % N3 == 0:
			time = int(100.*float(arrayz[0]))/100.
			print(' * t = '+str(time)+' fs')
			N = int(counter/N3)
			for l in range(0,N3):
				linex  = filex.readline()
				linex  = linex.strip()
				arrayx = linex.split()
				px.append(float(arrayx[3]))  				
			for i in range(0,N1):
				for k in range(0,N2):
					Z[i][k]= z[(N-1)*N3+i*N2+k]
					X[i][k]= x[(N-1)*N3+i*N2+k]
					p[i][k]= np.sqrt((pz[(N-1)*N3+i*N2+k]**2.)+(px[(N-1)*N3+i*N2+k]**2.))
			Maxval = np.matrix(p).max()
			Power=np.int(np.log(Maxval)/np.log(10.))
			P = p / (10.**Power)
			cmap = plt.get_cmap(cmap)
			norm = cm.colors.Normalize(vmax=np.matrix(P).max(), vmin=np.matrix(P).min())
			fig=plt.figure()
			plt.rc('text', usetex=True)
			plt.pcolormesh(Z,X,P,cmap=cmap,norm=norm,vmax=np.matrix(P).max(), vmin=np.matrix(P).min())
			cbar=plt.colorbar()
			cbar.ax.tick_params(labelsize=16) 
			plt.gca().set_aspect('equal')
			#font = {'family': 'non-serif',
			font = {
					'style':  'normal',
        			'color':  'black',
       		 		'weight': 'normal',
       		 		'size': 16,
       		 	}
			plt.title(str('%.0e' % (10**(-Power)))+r'$\times$'+title+' at '+str(time)+' fs', fontdict=font)
			plt.xticks(fontsize=16)
			plt.xlabel(r'$z\,(\mu\mathrm{m})$', fontdict=font)
			plt.xlim([np.amin(z[(N-1)*N3:N*N3]),np.amax(z[(N-1)*N3:N*N3])])
			plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=font)
			plt.yticks(fontsize=16)
			plt.ylim([np.amin(x[(N-1)*N3:N*N3]),np.amax(x[(N-1)*N3:N*N3])])
			fig.savefig(name+str(N)+'.png',bbox_inches='tight')
			plt.close(fig)
	filez.close()
	filex.close()
	return;
	
def read_and_plot_distribution(Ne,e0,mu,Nmu,filename_psi0,filename_psi1x,filename_psi1z,name):
	Ne     = int(Ne)
	Nmu    = int(Nmu)
	NeNmu  = Ne*Nmu
	Ntheta = 500
	theta  = np.linspace(0., 2. * np.pi, Ntheta)
	g  = np.zeros(Ne)
	p  = np.zeros(Ne)
	v  = np.zeros(Ne)
	Pz = np.zeros((Ntheta,Ne))
	Px = np.zeros((Ntheta,Ne))
	c   = 2.9979e8
	m   = 9.1091e-31
	keV = 1.6022e-16
	mc2 = (m*c)**2.
	mc3 = (m*c)**3.
	unit = (1./keV)*(c/mc2)*mc3
	for l in range(0,Ne):
		g[l]=1. + (e0[l]/511.)
		p[l]=np.sqrt((g[l]**2.)-1.)
		v[l]=np.sqrt(1.-(g[l]**(-2.)))
		for t in range(0,Ntheta):
			Pz[t][l]=p[l]*math.cos(theta[t])
			Px[t][l]=p[l]*math.sin(theta[t])
	time   = []
	z      = []
	x      = []
	f0     = []
	f1x    = []
	f1z    = []
	Ox = np.zeros(Ne)
	Oz = np.zeros(Ne)
	O  = np.zeros(Ne)
	ax = np.zeros(Ne)
	az = np.zeros(Ne)
	a  = np.zeros(Ne)
	norma = np.zeros(Ne)
	F = np.zeros((Ntheta,Ne))
	psi0_file  = open(filename_psi0,  'r')
	psi1x_file = open(filename_psi1x, 'r')
	psi1z_file = open(filename_psi1z, 'r')
	counter = 0
	if mu == 1: #z
		colx = 3
		colz = 1
	elif mu == 2: #x
		colx = 1
		colz = 3
	for line in psi0_file:
		line      = line.strip()
		array     = line.split()
		time.append(float(array[0]))
		x.append(float(array[colx]))
		z.append(float(array[colz]))
		f0.append(float(array[4]))
		counter = counter + 1
		if (counter % NeNmu == 0):
			N = int(counter/NeNmu)
			print(' * t = '+str(time[(N-1)*NeNmu])+' fs')
			for k in range(0,NeNmu):
				line_x  = psi1x_file.readline()
				line_x  = line_x.strip()
				array_x = line_x.split()
				f1x.append(float(array_x[4]))
				line_z  = psi1z_file.readline()
				line_z  = line_z.strip()
				array_z = line_z.split()
				f1z.append(float(array_z[4]))	
			#for i in range(1,Nmu+1):
			i =1 + int(Nmu/2.)
			for l in range(0,Ne):
				n = (N-1)*NeNmu+(i-1)*Ne+l
				if f0[n]!=0.:
					Ox[l]=f1x[n]/f0[n]
					Oz[l]=f1z[n]/f0[n]
				else:
					Ox[l]=0.
					Oz[l]=0.
				O[l]=np.sqrt((Ox[l]**2.)+(Oz[l]**2.))
				if O[l]>=0.998585:
					O[l]=0.998585
				ax[l]=3.*Ox[l]/(1.-(0.5*(O[l]**2.))*(1.+(O[l]**2.)))
				az[l]=3.*Oz[l]/(1.-(0.5*(O[l]**2.))*(1.+(O[l]**2.)))
				a[l]=np.sqrt((ax[l]**2.)+(az[l]**2.))
				if a[l]==0.:
					norma[l]=1./(4.*np.pi)
				elif a[l]>709.:
					norma[l]=1.373e-306
				else:
					norma[l]=a[l]/(4.*np.pi*np.sinh(a[l]))
				norma[l] = f0[n]*norma[l]*(v[l]/(p[l]**2.))*unit
				for t in range(0,Ntheta):
					if ax[l]*math.sin(theta[t])+az[l]*math.cos(theta[t])>709.:
						F[t][l]=norma[l]*8.2184e307
					else:
						F[t][l]=norma[l]*math.exp(ax[l]*math.sin(theta[t])+az[l]*math.cos(theta[t]))
					F[t][l]=np.log(1.+F[t][l])/np.log(10.)
			dir= os.path.dirname(name+str(i)+'/')
			if not os.path.exists(dir):
				os.mkdir(dir)
			if mu == 1: #z
				fb_mu = '/fb_z'
			elif mu ==2: #x
				fb_mu = '/fb_x'
			fig=plt.figure()
			plt.rc('text', usetex=True)
			cmap = plt.get_cmap('nipy_spectral')
			norm = cm.colors.Normalize(vmax=np.amax(F), vmin=np.amin(F))
			plt.pcolormesh(Pz,Px,F,cmap=cmap,norm=norm,vmax=np.amax(F), vmin=np.amin(F))
			plt.colorbar()
			font = {'family': 'serif',
					'style':  'normal',
					'weight': 'normal',
					'size': 18
					}
			plt.title(r'$\log_{10}(\mathbf{f_b}(x_0,\,z_0,\,\mathbf{p},\,t=$'
			+str(int(time[n]))+r'$\mathrm{fs})\,(\mathrm{cm}^{-3}.\mathrm{m_e c}^{-3}))$',
			fontdict=font)
			plt.xticks(fontsize=18)
			plt.xlabel(r'$p_z\,(m_e c)$', fontdict=font)
			plt.xlim([np.amin(Pz)-1,np.amax(Pz)+1])
			plt.ylabel(r'$p_x\,(m_e c)$', fontdict=font)
			plt.yticks(fontsize=18)
			plt.ylim([np.amin(Px)-1.,np.amax(Px)+1.])
			pos_x1 = np.amin(Pz)-0.5
			pos_y1 = np.amax(Px)*(10./11.)
			pos_x2 = 0.5*(np.amax(Pz)+0.5)
			pos_y2 = pos_y1
			plt.text(pos_x1, pos_y1, r'$x_0=$'+str(np.floor(100*x[n])/100)+r'$\mu\mathrm{m}$', color='red',fontsize=18)
			plt.text(pos_x2, pos_y2, r'$z_0=$'+str(np.floor(100*z[n])/100)+r'$\mu\mathrm{m}$', color='red',fontsize=18)
			fig.savefig(dir+fb_mu+str(i)+'_'+str(N)+'.png',bbox_inches='tight')
			plt.close(fig)
	psi0_file.close()
	psi1x_file.close()
	psi1z_file.close()
