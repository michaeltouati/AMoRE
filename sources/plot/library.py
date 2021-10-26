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
Python functions library for plotting AMoRE simulation results
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

FONT_SIZE = 16

FONT = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': FONT_SIZE,
        }

def get_results_dir():
    """
    Return AMoRE simulation name
    """
    with open('input-deck', 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            if array[0] == '#simu' :
                name = array[1]
                break
    to_print      = ' ' + name + ' AMoRE SIMULATION PLOTS :'
    line_to_print = ' '
    for _ in range(0,len(to_print)-1):
        line_to_print += '='
    print(str(line_to_print))
    print(to_print)
    print(str(line_to_print))
    return name

def create_dir(name):
    """
    Create directory entitled name
    """
    dir_name= os.path.dirname(name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def find_spatial_simulation_box_dimension(file_name):
    """
    Return AMoRE simulation spatial bins numbers
    """
    t_0  = []
    x_0  = []
    with open(file_name, 'r', encoding='utf-8') as nb_file0 :
        line = nb_file0.readline()
        line = line.strip()
        array = line.split()
        t_0.append(float(array[0]))
        x_0.append(float(array[1]))
        counter = 0
        for line in nb_file0:
            line      = line.strip()
            array     = line.split()
            t_0.append(float(array[0]))
            x_0.append(float(array[1]))
            counter = counter + 1
            if x_0[counter]!=x_0[counter-1]:
                n_z = counter
                break
        for line in nb_file0:
            line      = line.strip()
            array     = line.split()
            t_0.append(float(array[0]))
            counter = counter + 1
            if t_0[counter]!=t_0[counter-1]:
                n_x_n_z = counter
                break
    n_x = n_x_n_z / n_z
    return [int(n_z),int(n_x)]

def find_energy_bins(file_name):
    """
    Return AMoRE simulation beam electrons kinetic energy bins number
    """
    x_0  = []
    e_0  = []
    with open(file_name, 'r', encoding='utf-8') as psi0_file0 :
        line = psi0_file0.readline()
        line = line.strip()
        array = line.split()
        x_0.append(float(array[1]))
        e_0.append(float(array[2]))
        counter = 0
        for line in psi0_file0:
            line      = line.strip()
            array     = line.split()
            x_0.append(float(array[1]))
            e_0.append(float(array[2]))
            counter = counter + 1
            if x_0[counter]!=x_0[counter-1]:
                n_eps = counter
                break
    return [int(n_eps),np.array(e_0)]

def read_and_plot_curve(filename,xlabel,ylabel,title,name,logx,logy) :
    """
    Read and plot AMoRE scalar simulation result file entitled filename
    """
    [x_plt,y_plt] = read_file_and_define_two_first_cols(filename)
    fig=plt.figure()
    plt.rc('text', usetex=True)
    if logx==1:
        if logy ==0:
            plt.semilogx(x_plt, y_plt,linewidth=2)
        else:
            plt.loglog(x_plt, y_plt,linewidth=2)
    plt.title(title, fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(xlabel, fontdict=FONT)
    plt.ylabel(ylabel, fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    plt.grid(which='both', axis='both')
    fig.savefig(name,bbox_inches='tight')
    plt.close(fig)

def read_and_plot_two_log_curve(filename1,filename2,
                                label1,label2,
                                loc,xlabel,ylabel,title,
                                name):
    """
    Read and plot two AMoRE scalar simulation result files entitled filename1(2)
    """
    [x_1,y_1] = read_file_and_define_two_first_cols(filename1)
    [x_2,y_2] = read_file_and_define_two_first_cols(filename2)
    fig=plt.figure()
    plt.rc('text', usetex=True)
    plt.loglog(x_1, y_1,'red',linewidth=2,label=label1)
    plt.loglog(x_2, y_2,'blue',linewidth=2,label=label2)
    plt.title(title, fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(xlabel, fontdict=FONT)
    plt.ylabel(ylabel, fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    leg = plt.legend(loc=loc,fontsize=FONT_SIZE, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    plt.grid(which='both', axis='both')
    fig.savefig(name,bbox_inches='tight')
    plt.close(fig)

def read_file_and_define_first_col(filename):
    """
    Read and return the first column of an
    AMoRE simulation result file entitled filename
    """
    vector = []
    with open(filename, 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            vector.append(float(array[0]))
    return np.array(vector)

def read_file_and_define_second_col(filename):
    """
    Read and return the second column of an
    AMoRE simulation result file entitled filename
    """
    vector = []
    with open(filename, 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            vector.append(float(array[1]))
    return np.array(vector)

def read_file_and_define_two_first_cols(filename):
    """
    Read and return the two first columns of
    an AMoRE simulation result file entitled filename
    """
    x_0 = []
    y_0 = []
    with open(filename, 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            x_0.append(float(array[0]))
            y_0.append(float(array[1]))
    return [x_0,y_0]

def read_and_plot_2d_pcolormesh(filename,n_1,n_2,cmap,title,name,log):
    """
    Read and plot an AMoRE 2D map simulation result file entitled filename
    """
    n_1   = int(n_1)
    n_2   = int(n_2)
    n_3   = n_1*n_2
    x_plt = []
    z_plt = []
    p_plt = []
    z_map = np.zeros((n_1,n_2))
    x_map = np.zeros((n_1,n_2))
    p_map = np.zeros((n_1,n_2))
    with open(filename, 'r', encoding='utf-8') as file :
        counter = 0
        for line in file:
            line      = line.strip()
            array     = line.split()
            x_plt.append(float(array[1]))
            z_plt.append(float(array[2]))
            if log==1:
                p_plt.append(np.log(float(array[3]))/np.log(10))
            else:
                p_plt.append(float(array[3]))
            counter = counter + 1
            if counter % n_3 == 0:
                time     = float(array[0])
                str_time = f"{time:1.4E}"
                print(' * t = '+str_time+' fs')
                n_t   = int(counter / n_3)
                val_a = abs(np.amax(p_plt[(n_t-1)*n_3:n_t*n_3]))
                val_b = abs(np.amin(p_plt[(n_t-1)*n_3:n_t*n_3]))
                val   = max(val_a,val_b)
                if (val != 0.) and (log!=1):
                    power_log10=np.floor(np.log(val)/np.log(10.))
                else:
                    power_log10 = 0.
                for index in range((n_t-1)*n_3,n_t*n_3):
                    p_plt[index] = p_plt[index] / (10.**power_log10)
                for i in range(0,n_1):
                    for k in range(0,n_2):
                        z_map[i][k]=z_plt[(n_t-1)*n_3+i*n_2+k]
                        x_map[i][k]=x_plt[(n_t-1)*n_3+i*n_2+k]
                        p_map[i][k]=p_plt[(n_t-1)*n_3+i*n_2+k]
                if cmap == 'seismic' :
                    max_val_a = np.amax(p_plt[(n_t-1)*n_3:n_t*n_3])
                    max_val_b = abs(np.amin(p_plt[(n_t-1)*n_3:n_t*n_3]))
                    max_val = max(max_val_a,max_val_b)
                    min_val = -max_val
                else:
                    max_val = np.amax(p_plt[(n_t-1)*n_3:n_t*n_3])
                    min_val = np.amin(p_plt[(n_t-1)*n_3:n_t*n_3])
                norm = cm.colors.Normalize(vmax=max_val, vmin=min_val)
                fig=plt.figure()
                plt.rc('text', usetex=True)
                plt.pcolormesh(z_map,x_map,p_map,
                    cmap=plt.get_cmap(cmap),norm=norm,shading='gouraud')
                cbar=plt.colorbar()
                cbar.ax.tick_params(labelsize=FONT_SIZE)
                plt.gca().set_aspect('equal')
                if log!=1:
                    plot_title  = f"{10**(-power_log10):.0E}"+r'$\times$'
                    plot_title += title+'\n at '+str_time+' fs'
                    plt.title(plot_title, fontdict=FONT)
                else:
                    plt.title(title+'\n at '+str_time+' fs', fontdict=FONT)
                plt.xticks(fontsize=FONT_SIZE)
                plt.xlabel(r'$z\,(\mu\mathrm{m})$', fontdict=FONT)
                plt.xlim([np.amin(z_plt[(n_t-1)*n_3:n_t*n_3]),np.amax(z_plt[(n_t-1)*n_3:n_t*n_3])])
                plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=FONT)
                plt.yticks(fontsize=FONT_SIZE)
                plt.ylim([np.amin(x_plt[(n_t-1)*n_3:n_t*n_3]),np.amax(x_plt[(n_t-1)*n_3:n_t*n_3])])
                fig.savefig(name+str(n_t)+'.png',bbox_inches='tight')
                plt.close(fig)

def read_and_plot_2d_pcolormesh_abs(filenamez,filenamex,n_1,n_2,cmap,title,name):
    """
    Read and plot the absolute value of a vactor from its two
    AMoRE 2D map simulation result file entitled filename1(2)
    """
    n_1 = int(n_1)
    n_2 = int(n_2)
    n_3 = n_1*n_2
    x_plt = []
    z_plt = []
    p_z = []
    p_x = []
    z_map = np.zeros((n_1,n_2))
    x_map = np.zeros((n_1,n_2))
    p_map = np.zeros((n_1,n_2))
    with open(filenamez, 'r', encoding='utf-8') as filez :
        with open(filenamex, 'r', encoding='utf-8') as filex :
            counter = 0
            for linez in filez:
                linez      = linez.strip()
                arrayz     = linez.split()
                x_plt.append(float(arrayz[1]))
                z_plt.append(float(arrayz[2]))
                p_z.append(float(arrayz[3]))
                counter = counter + 1
                if counter % n_3 == 0:
                    time     = float(arrayz[0])
                    str_time = f"{time:1.4E}"
                    print(' * t = '+str_time+' fs')
                    n_t = int(counter/n_3)
                    for _ in range(0,n_3):
                        linex  = filex.readline()
                        linex  = linex.strip()
                        arrayx = linex.split()
                        p_x.append(float(arrayx[3]))
                    for i in range(0,n_1):
                        for k in range(0,n_2):
                            z_map[i][k] = z_plt[(n_t-1)*n_3+i*n_2+k]
                            x_map[i][k] = x_plt[(n_t-1)*n_3+i*n_2+k]
                            p_map_z2    = p_z[(n_t-1)*n_3+i*n_2+k]**2.
                            p_map_x2    = p_x[(n_t-1)*n_3+i*n_2+k]**2.
                            p_map[i][k] = np.sqrt(p_map_z2+p_map_x2)
                    max_val = np.matrix(p_map).max()
                    power_log10=np.int(np.log(max_val)/np.log(10.))
                    p_map = p_map / (10.**power_log10)
                    cmap = plt.get_cmap(cmap)
                    max_val = np.matrix(p_map).max()
                    min_val = np.matrix(p_map).min()
                    norm = cm.colors.Normalize(vmax=max_val, vmin=min_val)
                    fig=plt.figure()
                    plt.rc('text', usetex=True)
                    plt.pcolormesh(z_map,x_map,p_map,cmap=cmap,norm=norm,shading='gouraud')
                    cbar=plt.colorbar()
                    cbar.ax.tick_params(labelsize=FONT_SIZE)
                    plt.gca().set_aspect('equal')
                    plot_title  = f"{10**(-power_log10):.0E}"
                    plot_title += r'$\times$'+title+' at '+str_time+' fs'
                    plt.title(plot_title, fontdict=FONT)
                    plt.xticks(fontsize=FONT_SIZE)
                    plt.xlabel(r'$z\,(\mu\mathrm{m})$', fontdict=FONT)
                    z_plt_min = np.amin(z_plt[(n_t-1)*n_3:n_t*n_3])
                    z_plt_max = np.amax(z_plt[(n_t-1)*n_3:n_t*n_3])
                    plt.xlim([z_plt_min,z_plt_max])
                    plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=FONT)
                    plt.yticks(fontsize=FONT_SIZE)
                    x_plt_min = np.amin(x_plt[(n_t-1)*n_3:n_t*n_3])
                    x_plt_max = np.amax(x_plt[(n_t-1)*n_3:n_t*n_3])
                    plt.ylim([x_plt_min,x_plt_max])
                    fig.savefig(name+str(n_t)+'.png',bbox_inches='tight')
                    plt.close(fig)

def read_and_plot_distribution(**kwargs):
    """
    Plot the beam electron distribution from the
    AMoRE elecron beam distribution angular moment simulation result files entitled
    filename_psi0, filename_psi1x, filename_psi1z
    """
    n_e_n_mu_grk  = kwargs['n_e']*kwargs['n_mu_grk']
    n_theta       = 500
    theta         = np.linspace(0., 2. * np.pi, n_theta)
    gamma_l       = np.zeros(kwargs['n_e'])
    p_b           = np.zeros(kwargs['n_e'])
    v_b           = np.zeros(kwargs['n_e'])
    p_z           = np.zeros((n_theta,kwargs['n_e']))
    p_x           = np.zeros((n_theta,kwargs['n_e']))
    c_light       = 2.9979e8
    m_ele         = 9.1091e-31
    kev_in_joules = 1.6022e-16
    mc2           = (m_ele*c_light)**2.
    mc3           = (m_ele*c_light)**3.
    unit = (1./kev_in_joules)*(c_light/mc2)*mc3
    for l_eps in range(0,kwargs['n_e']):
        gamma_l[l_eps]=1. + (kwargs['e_0'][l_eps]/511.)
        p_b[l_eps]=np.sqrt((gamma_l[l_eps]**2.)-1.)
        v_b[l_eps]=np.sqrt(1.-(gamma_l[l_eps]**(-2.)))
        for i_theta in range(0,n_theta):
            p_z[i_theta][l_eps]=p_b[l_eps]*np.cos(theta[i_theta])
            p_x[i_theta][l_eps]=p_b[l_eps]*np.sin(theta[i_theta])
    time  = []
    z_0   = []
    x_0   = []
    f_0   = []
    f1x   = []
    f1z   = []
    o_x   = np.zeros(kwargs['n_e'])
    o_z   = np.zeros(kwargs['n_e'])
    o_m   = np.zeros(kwargs['n_e'])
    a_x   = np.zeros(kwargs['n_e'])
    a_z   = np.zeros(kwargs['n_e'])
    a_m   = np.zeros(kwargs['n_e'])
    norma = np.zeros(kwargs['n_e'])
    distr = np.zeros((n_theta,kwargs['n_e']))
    with open(kwargs['filename_psi0'],  'r', encoding='utf-8') as psi0_file :
        with open(kwargs['filename_psi1x'], 'r', encoding='utf-8') as psi1x_file :
            with open(kwargs['filename_psi1z'], 'r', encoding='utf-8') as psi1z_file :
                counter = 0
                if kwargs['mu_grk'] == 1: #z
                    colx = 3
                    colz = 1
                    fb_mu = '/fb_z'
                elif kwargs['mu_grk'] == 2: #x
                    colx = 1
                    colz = 3
                    fb_mu = '/fb_x'
                for line in psi0_file:
                    line      = line.strip()
                    array     = line.split()
                    time.append(float(array[0]))
                    x_0.append(float(array[colx]))
                    z_0.append(float(array[colz]))
                    f_0.append(float(array[4]))
                    counter = counter + 1
                    if counter % n_e_n_mu_grk == 0 :
                        n_t = int(counter/n_e_n_mu_grk)
                        str_time = f"{time[(n_t-1)*n_e_n_mu_grk]:1.4E}"
                        print(' * t = '+str_time+' fs')
                        for _ in range(0,n_e_n_mu_grk):
                            line_x  = psi1x_file.readline()
                            line_x  = line_x.strip()
                            array_x = line_x.split()
                            f1x.append(float(array_x[4]))
                            line_z  = psi1z_file.readline()
                            line_z  = line_z.strip()
                            array_z = line_z.split()
                            f1z.append(float(array_z[4]))
                        i =1 + int(kwargs['n_mu_grk']/2.)
                        for l_eps in range(0,kwargs['n_e']):
                            n_0 = (n_t-1)*n_e_n_mu_grk+(i-1)*kwargs['n_e']+l_eps
                            if f_0[n_0]!=0.:
                                o_x[l_eps]=f1x[n_0]/f_0[n_0]
                                o_z[l_eps]=f1z[n_0]/f_0[n_0]
                            else:
                                o_x[l_eps]=0.
                                o_z[l_eps]=0.
                            o_m[l_eps]=np.sqrt((o_x[l_eps]**2.)+(o_z[l_eps]**2.))
                            if o_m[l_eps]>=0.998585:
                                o_m[l_eps]=0.998585
                            o_m_2 = o_m[l_eps]**2.
                            a_x[l_eps]=3.*o_x[l_eps]/(1.-(0.5*o_m_2)*(1.+o_m_2))
                            a_z[l_eps]=3.*o_z[l_eps]/(1.-(0.5*o_m_2)*(1.+o_m_2))
                            a_m[l_eps]=np.sqrt((a_x[l_eps]**2.)+(a_z[l_eps]**2.))
                            if a_m[l_eps]==0.:
                                norma[l_eps]=1./(4.*np.pi)
                            elif a_m[l_eps]>709.:
                                norma[l_eps]=1.373e-306
                            else:
                                norma[l_eps]=a_m[l_eps]/(4.*np.pi*np.sinh(a_m[l_eps]))
                            norma[l_eps] = f_0[n_0]*norma[l_eps]*(v_b[l_eps]/(p_b[l_eps]**2.))*unit
                            for i_theta in range(0,n_theta):
                                a_thres  = a_x[l_eps]*np.sin(theta[i_theta])
                                a_thres += a_z[l_eps]*np.cos(theta[i_theta])
                                if a_thres>709.:
                                    distr[i_theta][l_eps]=norma[l_eps]*8.2184e307
                                else:
                                    phi_var  = a_x[l_eps]*np.sin(theta[i_theta])
                                    phi_var += a_z[l_eps]*np.cos(theta[i_theta])
                                    distr[i_theta][l_eps] = norma[l_eps]*np.exp(phi_var)
                                distr[i_theta][l_eps]=np.log(1.+distr[i_theta][l_eps])/np.log(10.)
                        create_dir(kwargs['name']+str(i)+'/')
                        fig=plt.figure()
                        plt.rc('text', usetex=True)
                        cmap = plt.get_cmap('nipy_spectral')
                        distr_min = np.matrix(distr).min()
                        distr_max = np.matrix(distr).max()
                        norm = cm.colors.Normalize(vmax=distr_max, vmin=distr_min)
                        plt.pcolormesh(p_z,p_x,distr,cmap=cmap,norm=norm,shading='gouraud')
                        plt.colorbar()
                        plot_title  = r'$\log_{10}(f_b(x_0,\,z_0,\,p_x,\,p_{y_0},\,p_z,\,t=$'
                        plot_title += f"{time[n_0]:1.4E}"
                        plot_title += r'$\mathrm{fs})\,$'
                        plot_title += r'$(\mathrm{cm}^{-3}.\mathrm{m_e c}^{-3}))$'
                        plt.title(plot_title,fontdict=FONT)
                        plt.xticks(fontsize=FONT_SIZE)
                        p_z_min = np.matrix(p_z).min()
                        p_z_max = np.matrix(p_z).max()
                        plt.xlabel(r'$p_z\,(m_e c)$', fontdict=FONT)
                        plt.xlim([p_z_min-1,p_z_max+1])
                        plt.ylabel(r'$p_x\,(m_e c)$', fontdict=FONT)
                        plt.yticks(fontsize=FONT_SIZE)
                        p_x_min = np.matrix(p_x).min()
                        p_x_max = np.matrix(p_x).max()
                        plt.ylim([p_x_min-1.,p_x_max+1.])
                        pos_x = p_z_min-0.5
                        pos_y = p_x_min-1. + 0.875*(p_x_max-p_x_min+2.)
                        txt_plt  =  r'$x_0\,=$'
                        txt_plt +=str(np.floor(100*x_0[n_0])/100)+r'$\,\mu\mathrm{m},\,$'
                        txt_plt += r'$z_0=$'+str(np.floor(100*z_0[n_0])/100)+r'$\,\mu\mathrm{m}$'
                        txt_plt += ' and \n'+r'$p_{y_0}=0\,\mathrm{m_e c}$'
                        plt.text(pos_x, pos_y, txt_plt,  color='red',fontsize=FONT_SIZE)
                        file_name = kwargs['name']+str(i)+'/'+fb_mu+str(i)+'_'+str(n_t)+'.png'
                        fig.savefig(file_name,bbox_inches='tight')
                        plt.close(fig)
