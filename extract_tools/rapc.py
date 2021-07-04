import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import os

def read_and_plot_curve(filename,xlabel,ylabel,title,name,logx,logy):
	x = []
	y = []
	file = open(filename, 'r')
	for line in file:
		line      = line.strip()
		array     = line.split()
		x.append(float(array[0]))
		y.append(float(array[1]))
	file.close()   
	fig=plt.figure()
	plt.rc('text', usetex=True)
	if logx==1:
		if logy ==0:
			plt.semilogx(x, y,linewidth=2)	
		else:
			plt.loglog(x, y,linewidth=2)
	#font = {'family': 'non-serif',
	#		'style':  'normal',
	font = {'style':  'normal',
        	'color':  'black',
    	    'weight': 'normal',
        	'size': 16,
       	 }
	plt.title(title, fontdict=font)
	plt.xticks(fontsize=16)
	plt.xlabel(xlabel, fontdict=font)
	plt.ylabel(ylabel, fontdict=font)
	plt.yticks(fontsize=16)
	fig.savefig(name,bbox_inches='tight')
	plt.close(fig)
	return;

def read_and_plot_two_log_curve(filename1,filename2,label1,label2,loc,xlabel,ylabel,title,name):
	x1 = []
	y1 = []
	file1 = open(filename1, 'r')
	for line in file1:
		line      = line.strip()
		array     = line.split()
		x1.append(float(array[0]))
		y1.append(float(array[1]))
	file1.close() 
	x2 = []
	y2 = []
	file2 = open(filename2, 'r')
	for line in file2:
		line      = line.strip()
		array     = line.split()
		x2.append(float(array[0]))
		y2.append(float(array[1]))
	file2.close()   
	fig=plt.figure()
	plt.rc('text', usetex=True)
	plt.loglog(x1, y1,'red',linewidth=2,label=label1)
	plt.loglog(x2, y2,'blue',linewidth=2,label=label2)
	#font = {'family': 'non-serif',
	font = {'style':  'normal',
        	'color':  'black',
    	    'weight': 'normal',
        	'size': 16,
       	 }
	plt.title(title, fontdict=font)
	plt.xticks(fontsize=16)
	plt.xlabel(xlabel, fontdict=font)
	plt.ylabel(ylabel, fontdict=font)
	plt.yticks(fontsize=16)
	leg = plt.legend(loc=loc,fontsize=16, fancybox=True)
	leg.get_frame().set_alpha(0.5)
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
			time = int(float(array[0]))
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
			time = int(float(arrayz[0]))
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
			Power=math.floor(np.log(Maxval)/np.log(10.))
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
			print('-at t = '+str(time[(N-1)*NeNmu])+' fs')
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
