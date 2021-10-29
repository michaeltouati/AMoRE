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

def get_log_axis_min_value(vector):
    """
    Return the maximum value 10^n with n integer
    lower than all vector values
    """
    vector_min = np.log(np.min(vector))/np.log(10.)
    vector_min = 10.**np.floor(vector_min)
    return vector_min

def get_log_axis_max_value(vector):
    """
    Return the minimum value 10^n with n integer
    greater than all vector values
    """
    vector_max = np.log(np.max(vector))/np.log(10.)
    vector_max = 10.**np.ceil(vector_max)
    return vector_max

def make_scalars_plot_figure(**kwargs):
    """
    Plot and save a .png image of scalar fields (maximum 10)
    kwargs keys :
    * scalars  : xplot_min, xplot_max,
                 yplot_min, yplot_max
    * vectors  : xplot,
                 yplot1, ... yplot10
    * strings  : xlabel, ylabel,
                 color1, ..., color10,
                 legend1, ..., legend9
                 title, filename
    * logicals : logx, logy, grid
    """
    fig=plt.figure()
    plt.rc('text', usetex=True)
    for i in range(1,11):
        condition = ('yplot'+str(i) in kwargs)
        condition = condition and ('legend'+str(i) in kwargs)
        condition = condition and ('color'+str(i) in kwargs)
        if condition :
            yplot  = kwargs['yplot' +str(i)]
            legend = kwargs['legend'+str(i)]
            color  = kwargs['color' +str(i)]
            if kwargs['logx'] :
                if kwargs['logy'] :
                    plt.loglog(kwargs['xplot'],yplot,
                               linewidth=2,label=legend,color=color)
                else:
                    plt.semilogx(kwargs['xplot'],yplot,
                                 linewidth=2,label=legend,color=color)
            else :
                plt.plot(kwargs['xplot'],yplot,
                         linewidth=2,label=legend,color=color)
    leg = plt.legend(fontsize=FONT_SIZE,
                     fancybox=True,
                     bbox_to_anchor=[1., 1.],
                     loc='upper left')
    leg.get_frame().set_alpha(0.5)
    if ('xplot_min' in kwargs) and ('xplot_max' in kwargs) :
        plt.xlim([kwargs['xplot_min'],kwargs['xplot_max']])
    if ('yplot_min' in kwargs) and ('yplot_max' in kwargs) :
        plt.ylim([kwargs['yplot_min'],kwargs['yplot_max']])
    plt.title(kwargs['title'], fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(kwargs['xlabel'], fontdict=FONT)
    plt.ylabel(kwargs['ylabel'], fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    if kwargs['grid']:
        plt.grid(which='both', axis='both')
    fig.savefig(kwargs['filename'],bbox_inches='tight')
    plt.close(fig)

def make_scalar_plot_figure(**kwargs):
    """
    Plot and save a .png image of scalar fields (maximum 10)
    kwargs keys :
    * scalars  : xplot_min, xplot_max,
                 yplot_min, yplot_max
    * vectors  : xplot, yplot
    * strings  : xlabel, ylabel,
                 color
                 title, filename
    * logicals : logx, logy, grid
    """
    fig=plt.figure()
    plt.rc('text', usetex=True)
    if kwargs['logx'] :
        if kwargs['logy'] :
            plt.loglog(kwargs['xplot'],kwargs['yplot'],
                       linewidth=2,color=kwargs['color'])
        else :
            plt.semilogx(kwargs['xplot'],kwargs['yplot'],
                         linewidth=2,color=kwargs['color'])
    else :
        plt.plot(kwargs['xplot'],kwargs['yplot'],
                 linewidth=2,color=kwargs['color'])
    if ('xplot_min' in kwargs) and ('xplot_max' in kwargs) :
        plt.xlim([kwargs['xplot_min'],kwargs['xplot_max']])
    if ('yplot_min' in kwargs) and ('yplot_max' in kwargs) :
        plt.ylim([kwargs['yplot_min'],kwargs['yplot_max']])
    plt.title(kwargs['title'], fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(kwargs['xlabel'], fontdict=FONT)
    plt.ylabel(kwargs['ylabel'], fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    if kwargs['grid']:
        plt.grid(which='both', axis='both')
    fig.savefig(kwargs['filename'],bbox_inches='tight')
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
    return [np.array(x_0),np.array(y_0)]

def make_2d_field_pcolormesh_figure(**kwargs):
    """
    Plot and save a .png image of a 2D field map
    kwargs keys :
    * scalars   : xmap_min, xmap_max,
                  ymap_min, ymap_max,
                  zmap_min, zmap_max
    * 2D arrays : xmap, ymap, zmap
    * strings   : xlabel, ylabel,
                  colormap, title,
                  fig_file
    * logical   : eq_axis
    """
    cmap = plt.get_cmap(kwargs['colormap'])
    norm = cm.colors.Normalize(vmin=kwargs['zmap_min'],
                               vmax=kwargs['zmap_max'])
    fig=plt.figure()
    plt.rc('text', usetex=True)
    plt.pcolormesh(kwargs['xmap'],kwargs['ymap'],kwargs['zmap'],
                   cmap=cmap,norm=norm,shading='gouraud')
    cbar=plt.colorbar(format='%.2f')
    cbar.ax.tick_params(labelsize=FONT_SIZE)
    plt.title(kwargs['title'], fontdict=FONT)
    plt.xticks(fontsize=FONT_SIZE)
    plt.xlabel(kwargs['xlabel'], fontdict=FONT)
    plt.xlim([kwargs['xmap_min'],kwargs['xmap_max']])
    plt.ylabel(kwargs['ylabel'], fontdict=FONT)
    plt.yticks(fontsize=FONT_SIZE)
    plt.ylim([kwargs['ymap_min'],kwargs['ymap_max']])
    if kwargs['eq_axis'] :
        plt.gca().set_aspect('equal')
    fig.savefig(kwargs['fig_file'],bbox_inches='tight')
    plt.close(fig)

def mesh_grid(**kwargs):
    """
    Return 2D arrays from hydro./field AMoRE simulation results files
    kwargs keys :
    * 1D arrays : z_in, x_in, p_in
    * Integers  : n_0 (for t), n_1 (for x), n_2 (for z)
    """
    n_3   = kwargs['n_1']*kwargs['n_2']
    z_out = np.zeros((kwargs['n_1'],kwargs['n_2']))
    x_out = np.zeros((kwargs['n_1'],kwargs['n_2']))
    p_out = np.zeros((kwargs['n_1'],kwargs['n_2']))
    for i in range(0,kwargs['n_1']):
        for k in range(0,kwargs['n_2']):
            index       = (kwargs['n_0']-1)*n_3
            index      += i*kwargs['n_2']
            index      += k
            z_out[i][k] = kwargs['z_in'][index]
            x_out[i][k] = kwargs['x_in'][index]
            p_out[i][k] = kwargs['p_in'][index]
    return [z_out, x_out, p_out]

def normalize(**kwargs):
    """
    Normalize 2D array to plot and modify the title in consequence
    kwargs keys :
    * 1D array : p_in
    * Strings  : title, str_time
    * Logicals : log
    * Integers : n_0 (for t), n_1 (for x), n_2 (for z)
    """
    n_3       = kwargs['n_1']*kwargs['n_2']
    index_min = (kwargs['n_0']-1) * n_3
    index_max =  kwargs['n_0']    * n_3
    value     = max(abs(np.amax(kwargs['p_in'][index_min:index_max])),
                    abs(np.amin(kwargs['p_in'][index_min:index_max])))
    power_log10 = 0.
    plot_title  = ''
    if kwargs['log'] :
        plot_title  = kwargs['title']+'\n at '
        plot_title += kwargs['str_time']+' fs'
    else :
        if value != 0. :
            power_log10=np.floor(np.log(value)/np.log(10.))
        plot_title  = f"{10**(-power_log10):.0E}"+r'$\times$'
        plot_title += kwargs['title']+'\n at '
        plot_title += kwargs['str_time']+' fs'
    p_out = kwargs['p_in']
    for index in range(index_min,index_max):
        p_out[index] = kwargs['p_in'][index] / (10.**power_log10)
    return [plot_title, p_out]

def read_and_plot_2d_pcolormesh(**kwargs):
    """
    Read and plot an AMoRE 2D map simulation result file entitled filename
    kwargs keys  :
    * integers   : n_1, n_2
    * strings    : res_file, colormap, title, fig_file
    * Logicals   : log, eq_axis
    """
    x_plt = []
    z_plt = []
    p_plt = []
    with open(kwargs['res_file'], 'r', encoding='utf-8') as file :
        counter = 0
        for line in file:
            array = line.strip().split()
            x_plt.append(float(array[1]))
            z_plt.append(float(array[2]))
            if kwargs['log'] :
                p_plt.append(np.log(float(array[3]))/np.log(10))
            else:
                p_plt.append(float(array[3]))
            counter = counter + 1
            if counter % (kwargs['n_1']*kwargs['n_2']) == 0:
                print(' * t = '+f"{float(array[0]):1.4E}"+' fs')
                n_t = int( counter / (kwargs['n_1']*kwargs['n_2']) )
                [plot_title, p_plt]   = normalize(title    = kwargs['title'],
                                                  log      = kwargs['log'],
                                                  str_time = f"{float(array[0]):1.4E}",
                                                  p_in     = p_plt,
                                                  n_0      = n_t,
                                                  n_1      = kwargs['n_1'],
                                                  n_2      = kwargs['n_2'])
                [z_map, x_map, p_map] = mesh_grid(z_in = z_plt,
                                                  x_in = x_plt,
                                                  p_in = p_plt,
                                                  n_0  = n_t,
                                                  n_1  = kwargs['n_1'],
                                                  n_2  = kwargs['n_2'])
                if kwargs['colormap'] == 'seismic' :
                    max_val = max(np.matrix(p_map).max(),
                                  abs(np.matrix(p_map).min()))
                    min_val = -max_val
                else:
                    max_val = np.matrix(p_map).max()
                    min_val = np.matrix(p_map).min()
                make_2d_field_pcolormesh_figure(xmap     = z_map,
                                                ymap     = x_map,
                                                zmap     = p_map,
                                                colormap = kwargs['colormap'],
                                                xmap_min = np.matrix(z_map).min(),
                                                xmap_max = np.matrix(z_map).max(),
                                                ymap_min = np.matrix(x_map).min(),
                                                ymap_max = np.matrix(x_map).max(),
                                                zmap_min = min_val,
                                                zmap_max = max_val,
                                                xlabel   = kwargs['xlabel'],
                                                ylabel   = kwargs['ylabel'],
                                                title    = plot_title,
                                                eq_axis  = kwargs['eq_axis'],
                                                fig_file = kwargs['fig_file']+str(n_t)+'.png')

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
