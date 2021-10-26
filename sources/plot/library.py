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

font_size = 16

font = {'style':  'normal',
        'color':  'black',
        'weight': 'normal',
        'size': font_size,
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
    dir_name= os.path.dirname(name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def find_spatial_simulation_box_dimension(file_name):
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
    return [n_z,n_x]

def find_energy_bins(file_name):
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
    return [n_eps,np.array(e_0)]

def read_and_plot_curve(filename,xlabel,ylabel,title,name,logx,logy) :
    [x_plt,y_plt] = read_file_and_define_two_first_cols(filename)
    fig=plt.figure()
    plt.rc('text', usetex=True)
    if logx==1:
        if logy ==0:
            plt.semilogx(x_plt, y_plt,linewidth=2)
        else:
            plt.loglog(x_plt, y_plt,linewidth=2)
    plt.title(title, fontdict=font)
    plt.xticks(fontsize=font_size)
    plt.xlabel(xlabel, fontdict=font)
    plt.ylabel(ylabel, fontdict=font)
    plt.yticks(fontsize=font_size)
    plt.grid(which='both', axis='both')
    fig.savefig(name,bbox_inches='tight')
    plt.close(fig)

def read_and_plot_two_log_curve(filename1,filename2,
                                label1,label2,
                                loc,xlabel,ylabel,title,
                                name):
    [x_1,y_1] = read_file_and_define_two_first_cols(filename1)
    [x_2,y_2] = read_file_and_define_two_first_cols(filename2)
    fig=plt.figure()
    plt.rc('text', usetex=True)
    plt.loglog(x_1, y_1,'red',linewidth=2,label=label1)
    plt.loglog(x_2, y_2,'blue',linewidth=2,label=label2)
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

def read_file_and_define_first_col(filename):
    time = []
    with open(filename, 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            time.append(float(array[0]))
    return np.array(time)

def read_file_and_define_second_col(filename):
    vector = []
    with open(filename, 'r', encoding='utf-8') as file :
        for line in file:
            line      = line.strip()
            array     = line.split()
            vector.append(float(array[1]))
    return np.array(vector)

def read_file_and_define_two_first_cols(filename):
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
                cmap = plt.get_cmap(cmap)
                file_list = ('diag/jb_z[A_cm-2].dat',
                             'diag/je_z[A_cm-2].dat',
                             'diag/E_z[V_m-1].dat',
                             'diag/B_y[Tesla].dat')
                if filename in file_list :
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
                plt.pcolormesh(z_map,x_map,p_map,cmap=cmap,norm=norm,shading='gouraud')
                cbar=plt.colorbar()
                cbar.ax.tick_params(labelsize=font_size)
                plt.gca().set_aspect('equal')
                if log!=1:
                    plot_title = str('%.0e' % (10.**(-power_log10)))+r'$\times$'
                    plot_title = plot_title + title+' at '+str_time+' fs'
                    plt.title(plot_title, fontdict=font)
                else:
                    plt.title(title+' at '+str_time+' fs', fontdict=font)
                plt.xticks(fontsize=font_size)
                plt.xlabel(r'$z\,(\mu\mathrm{m})$', fontdict=font)
                plt.xlim([np.amin(z_plt[(n_t-1)*n_3:n_t*n_3]),np.amax(z_plt[(n_t-1)*n_3:n_t*n_3])])
                plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=font)
                plt.yticks(fontsize=font_size)
                plt.ylim([np.amin(x_plt[(n_t-1)*n_3:n_t*n_3]),np.amax(x_plt[(n_t-1)*n_3:n_t*n_3])])
                fig.savefig(name+str(n_t)+'.png',bbox_inches='tight')
                plt.close(fig)

def read_and_plot_2d_pcolormesh_abs(filenamez,filenamex,n_1,n_2,cmap,title,name):
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
                    cbar.ax.tick_params(labelsize=font_size)
                    plt.gca().set_aspect('equal')
                    plot_title = str('%.0e' % (10**(-power_log10)))
                    plot_title = plot_title +r'$\times$'+title+' at '+str_time+' fs'
                    plt.title(plot_title, fontdict=font)
                    plt.xticks(fontsize=font_size)
                    plt.xlabel(r'$z\,(\mu\mathrm{m})$', fontdict=font)
                    z_plt_min = np.amin(z_plt[(n_t-1)*n_3:n_t*n_3])
                    z_plt_max = np.amax(z_plt[(n_t-1)*n_3:n_t*n_3])
                    plt.xlim([z_plt_min,z_plt_max])
                    plt.ylabel(r'$x\,(\mu\mathrm{m})$', fontdict=font)
                    plt.yticks(fontsize=font_size)
                    x_plt_min = np.amin(x_plt[(n_t-1)*n_3:n_t*n_3])
                    x_plt_max = np.amax(x_plt[(n_t-1)*n_3:n_t*n_3])
                    plt.ylim([x_plt_min,x_plt_max])
                    fig.savefig(name+str(n_t)+'.png',bbox_inches='tight')
                    plt.close(fig)

def read_and_plot_distribution(n_e, e_0, mu_grk, n_mu_grk,
                               filename_psi0, filename_psi1x, filename_psi1z,
                               name):
    n_e           = int(n_e)
    n_mu_grk      = int(n_mu_grk)
    n_e_n_mu_grk  = n_e*n_mu_grk
    n_theta       = 500
    theta         = np.linspace(0., 2. * np.pi, n_theta)
    gamma_l       = np.zeros(n_e)
    p_b           = np.zeros(n_e)
    v_b           = np.zeros(n_e)
    p_z           = np.zeros((n_theta,n_e))
    p_x           = np.zeros((n_theta,n_e))
    c_light       = 2.9979e8
    m_ele         = 9.1091e-31
    kev_in_joules = 1.6022e-16
    mc2           = (m_ele*c_light)**2.
    mc3           = (m_ele*c_light)**3.
    unit = (1./kev_in_joules)*(c_light/mc2)*mc3
    for l_eps in range(0,n_e):
        gamma_l[l_eps]=1. + (e_0[l_eps]/511.)
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
    o_x   = np.zeros(n_e)
    o_z   = np.zeros(n_e)
    o_m   = np.zeros(n_e)
    a_x   = np.zeros(n_e)
    a_z   = np.zeros(n_e)
    a_m   = np.zeros(n_e)
    norma = np.zeros(n_e)
    distr = np.zeros((n_theta,n_e))
    with open(filename_psi0,  'r', encoding='utf-8') as psi0_file :
        with open(filename_psi1x, 'r', encoding='utf-8') as psi1x_file :
            with open(filename_psi1z, 'r', encoding='utf-8') as psi1z_file :
                counter = 0
                if mu_grk == 1: #z
                    colx = 3
                    colz = 1
                elif mu_grk == 2: #x
                    colx = 1
                    colz = 3
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
                        i =1 + int(n_mu_grk/2.)
                        for l_eps in range(0,n_e):
                            n_0 = (n_t-1)*n_e_n_mu_grk+(i-1)*n_e+l_eps
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
                        create_dir(name+str(i)+'/')
                        if mu_grk == 1: #z
                            fb_mu = '/fb_z'
                        elif mu_grk ==2: #x
                            fb_mu = '/fb_x'
                        fig=plt.figure()
                        plt.rc('text', usetex=True)
                        cmap = plt.get_cmap('nipy_spectral')
                        distr_min = np.matrix(distr).min()
                        distr_max = np.matrix(distr).max()
                        norm = cm.colors.Normalize(vmax=distr_max, vmin=distr_min)
                        plt.pcolormesh(p_z,p_x,distr,cmap=cmap,norm=norm,shading='gouraud')
                        plt.colorbar()
                        plot_title = r'$\log_{10}(\mathbf{f_b}(x_0,\,z_0,\,\mathbf{p},\,t=$'
                        plot_title = plot_title + str(int(time[n_0]))
                        plot_title = plot_title + r'$\mathrm{fs})\,$'
                        plot_title = plot_title + r'$(\mathrm{cm}^{-3}.\mathrm{m_e c}^{-3}))$'
                        plt.title(plot_title,fontdict=font)
                        plt.xticks(fontsize=font_size)
                        p_z_min = np.matrix(p_z).min()
                        p_z_max = np.matrix(p_z).max()
                        plt.xlabel(r'$p_z\,(m_e c)$', fontdict=font)
                        plt.xlim([p_z_min-1,p_z_max+1])
                        plt.ylabel(r'$p_x\,(m_e c)$', fontdict=font)
                        plt.yticks(fontsize=font_size)
                        p_x_min = np.matrix(p_x).min()
                        p_x_max = np.matrix(p_x).max()
                        plt.ylim([p_x_min-1.,p_x_max+1.])
                        pos_x_1 = p_z_min-0.5
                        pos_y_1 = p_x_max*(10./11.)
                        pos_x_2 = 0.5*(p_z_max+0.5)
                        pos_y_2 = pos_y_1
                        txt_xlabel = r'$x_0=$'+str(np.floor(100*x_0[n_0])/100)+r'$\mu\mathrm{m}$'
                        txt_ylabel = r'$z_0=$'+str(np.floor(100*z_0[n_0])/100)+r'$\mu\mathrm{m}$'
                        plt.text(pos_x_1, pos_y_1, txt_xlabel, color='red',fontsize=font_size)
                        plt.text(pos_x_2, pos_y_2, txt_ylabel, color='red',fontsize=font_size)
                        fig.savefig(name+str(i)+'/'+fb_mu+str(i)+'_'+str(n_t)+'.png',bbox_inches='tight')
                        plt.close(fig)
