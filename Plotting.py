# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 19:07:45 2016
Most recently edited: 16J27

@author: scott
"""

from matplotlib import pyplot as plt
import numpy as np
import os

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from EC import sync_metadata
else:                           #then we use relative import
    from .EC import sync_metadata



def plot_vs_potential(CV_and_MS, colors, tspan=0, RE_vs_RHE=None, A_el=None, 
                      ax1='new', ax2='new', overlay=0, logplot = [1,0], leg=1,
                      verbose=1):
    '''
    This will plot current and select MS signals vs E_we, as is the 
    convention for cyclic voltammagrams. added 16I29
    '''
    if verbose:
        print('\n\nfunction \'plot_vs_potential\' at your service!\n')
    
    #prepare axes
    if ax1 != 'new':
        figure1 = ax1.figure
    elif ax2 != 'new':
        figure1 = ax2.figure
    else:
        figure1 = plt.figure()
    if overlay:
        if ax1 == 'new':
            ax1 = figure1.add_subplot(111)
        if ax2 == 'new':
            ax2 = ax1.twinx()
    else:
        if ax1 == 'new':
            ax1 = figure1.add_subplot(211)
        if ax2 == 'new':
            ax2 = figure1.add_subplot(212)
    if type(logplot) is int:
        logplot = [logplot,logplot]
    if logplot[0]:
        ax1.set_yscale('log')
    if logplot[1]:
        ax2.set_yscale('log')
            
    # get EC data
    V_str, J_str = sync_metadata(CV_and_MS, RE_vs_RHE=RE_vs_RHE, A_el=A_el)
    V = CV_and_MS[V_str]
    J = CV_and_MS[J_str]

    #get time variable and plotting indexes
    t = CV_and_MS['time/s']
    if tspan == 0:                  #then use the whole range of overlap
        tspan = CV_and_MS['tspan_2']
    I_plot = np.array([i for (i,t_i) in enumerate(t) if tspan[0]<t_i and t_i<tspan[1]])
    
    #plot EC-lab data
    ax2.plot(V[I_plot],J[I_plot],'k-')      
        #maybe I should use EC.plot_cycles to have different cycles be different colors. Or rewrite that code here.
    ax2.set_xlabel(V_str)
    ax2.set_ylabel(J_str)
    
    for (mass, color) in colors.items():
        x_str = mass + '-x'
        y_str = mass + '-y'
        Y_str = mass + '-Y'     #-Y will be a QMS signal interpreted to the EC-lab time variable.
        x = CV_and_MS[x_str]
        y = CV_and_MS[y_str]
        Y = np.interp(t, x, y)  #obs! np.interp has a has a different argument order than Matlab's interp1
        CV_and_MS[Y_str] = Y    #add the interpolated value to the dictionary for future use
        ax1.plot(V[I_plot], Y[I_plot], color, label=mass)
    M_str = 'signal / [A]'
    ax1.set_xlabel(V_str)
    ax1.set_ylabel(M_str)
    if leg:
        ax1.legend()
    
    if verbose:
        print('\nfunction \'plot_vs_potential\' finished!\n\n')
        
        #parameter order of np.interp is different than Matlab's interp1
    return ax1, ax2, CV_and_MS    


def plot_vs_time(Dataset, cols_1='input', cols_2='input', verbose=1):
    '''
    Superceded by the more convenient plot_masses and plot_masses_and_I
    '''
    if verbose:
        print('\n\nfunction \'plot_vs_time\' at your service!')
    
    if cols_1=='input':
        data_cols = Dataset['data_cols']
        prompt = ('Choose combinations of time and non-time variables for axis 1, \n' +
            'with every other choice a time variable.')
        I_axis_1 = indeces_from_input(data_cols, prompt)
        cols_1 = [[data_cols[i], data_cols[j]] for i,j in zip(I_axis_1[::2],I_axis_1[1::2])]        
            
    figure1 = plt.figure()
    axes_1 = figure1.add_subplot(211)
    for pltpair in cols_1:
        label_object = pltpair[1][0:-2]
        if label_object:
            label_string = label_object.group()[:-1]
        else:
            label_string = pltpair[1]
        x = Dataset[pltpair[0]]
        y = np.log(Dataset[pltpair[1]])/np.log(10)
        axes_1.plot(x,y, label = label_string)
        
    axes_1.set_xlabel('time / s')
    axes_1.set_ylabel('log(signal/[a.u.])')
    axes_1.legend()    
    
    if cols_2=='input':
        
        data_cols = Dataset['data_cols']
        prompt = ('Choose combinations of time and non-time variables for axis 2, \n' +
            'with every other choice a time variable.')
        I_axis_2 = indeces_from_input(data_cols, prompt)
        cols_2 = [[data_cols[i], data_cols[j]] for i,j in zip(I_axis_2[::2],I_axis_2[1::2])]

    axes_2 = figure1.add_subplot(212)
    for pltpair in cols_2:
        label_string = pltpair[1]
        x = np.insert(Dataset[pltpair[0]],0,0)
        y = np.insert(Dataset[pltpair[1]],0,0)
        axes_2.plot(x,y,'k--',label=label_string)
    axes_2.set_ylabel('current / mA')
    axes_2.set_xlabel('time / s')
    axes_2.legend()
    #so capacitance doesn't blow it up:
    I_plt_top = np.where(x>2)[0][0]
    y_max = np.max(y[I_plt_top:])
    axes_2.set_ylim(np.min(y),y_max)
    if verbose:
        print('function \'plot_vs_time\' finished!\n\n')


def indeces_from_input(options, prompt):
    '''not sure I'll ever actually use this'''
    print(prompt + '\n... enter the indeces you\'re interested in, in order,' +
    'seperated by spaces, for example:\n>>>1 4 3')
    for nc, option in enumerate(options):
        print(str(nc) + '\t\t ' + options[nc])
    choice_string = input('\n')
    choices = choice_string.split(' ')
    choices = [int(choice) for choice in choices]
    return choices


def plot_masses(MS_Data, tspan=0, logplot=1, verbose=1,
                colors = {'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'}, 
                ax1='new', saveit=0, leg=1):
    '''
    plots selected masses for a selected time range from MS data or EC_MS data
    '''


    if verbose:
        print('\n\nfunction \'plot_masses\' at your service! \n Plotting from: ' + MS_Data['title'])

    if ax1 == 'new':
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)    
    lines = {}
    
    if tspan == 0:                  #then use the range of overlap
        tspan = CV_and_MS['tspan_2']    
    
    for mass, color in colors.items():
        if verbose:
            print('plotting: ' + mass)
        x = MS_Data[mass+'-x']
        y = MS_Data[mass+'-y']
        try:
            #np.where() keeps giving me headaches, so I'll try with list comprehension.
            index_list = np.array([i for (i,x_i) in enumerate(x) if tspan[0]<x_i and x_i<tspan[1]])     
            I_start = index_list[0]
            I_finish = index_list[-1]
            x = x[I_start:I_finish]
            y = y[I_start:I_finish]
        except IndexError:
            print('your tspan is probably fucked.\n x for ' + mass + ' goes from ' + str(x[0]) + ' to ' + str(x[-1]) +
                '\nand yet you ask for a tspan of ' + str(tspan[0]) + ' to ' + str(tspan[-1]))
        lines[mass] = ax1.plot(x, y, color, label = mass) 
        #as it is, lines is not actually used for anything         
    if leg:
        ax1.legend(loc = 'lower right')
    ax1.set_xlabel('time / [s]')
    y_string = 'signal / [A]'
    ax1.set_ylabel(y_string)           
    if logplot: 
        ax1.set_yscale('log') 
    if verbose:
        print('function \'plot_masses\' finsihed! \n\n')
    return ax1

def plot_masses_and_I(EC_and_MS, tspan=0, overlay=0, logplot=[1,0], verbose=1, 
                      colors={'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'}, 
                      plotpotential=1, RE_vs_RHE=None, A_el=None, 
                      saveit=0, title='default', leg=1):
    '''
    this plots current and potential on one axis and masses on another
    '''
    
    if verbose:
        print('\n\nfunction \'plot_masses_and_I\' at your service!\n Plotting from: ' + EC_and_MS['title'])
    
    figure1 = plt.figure()
    if overlay:
        ax1 = figure1.add_subplot(111)
        ax2 = ax1.twinx()
    else:
        ax1 = figure1.add_subplot(211)
        ax2 = figure1.add_subplot(212)
        
    if tspan == 0:                  #then use the whole range of overlap
        tspan = EC_and_MS['tspan_2']    
    plot_masses(EC_and_MS, tspan, logplot[0], verbose=verbose, colors=colors, ax1=ax1, saveit=0, leg=leg)
    
    t = EC_and_MS['time/s']
    
    V_str, J_str = sync_metadata(EC_and_MS, RE_vs_RHE=RE_vs_RHE, A_el=A_el) #added 16J27
        
    V = EC_and_MS[V_str]
    J = EC_and_MS[J_str]       
    
    ax2.plot(t, J, 'r')
    ax2.set_ylabel(J_str)
    ax2.set_xlabel('time / [s]')
    xlim = ax1.get_xlim()
    ax2.set_xlim(xlim)
    if logplot[1]: 
        ax2.set_yscale('log')  
    
    if plotpotential:
        ax3 = ax2.twinx()

        ax3.plot(t, V, 'k')
        ax3.set_ylabel(V_str)
        if len(logplot) >2:
            if logplot[2]:
                ax3.set_yscale('log')
        ax3.set_xlim(xlim)
    if saveit:
        if title == 'default':
            title == EC_and_MS['title'] + '.png'
        figure1.savefig(title)
        
    if verbose:
        print('function \'plot_masses_and_I\' finished!\n\n')
    if plotpotential:
        return ax1, ax2, ax3
    return ax1, ax2
    
    
    
if __name__ == '__main__':
    import os
    from Data_Importing import import_data
    from Combining import synchronize
    from EC import select_cycles, remove_delay
    
    plt.close('all')
    
    importrawdata = 1
    if importrawdata:
        default_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) 
        CV_data_0 = import_data(default_directory + os.sep) # + '18_CO_dose_and_strip_C01.mpt', data_type='EC')
        MS_data_0 = import_data(default_directory + os.sep,# + 'QMS_16I27_18h35m30.txt'
                                data_type='MS')
    
    CV_data = select_cycles(CV_data_0,[1,2])    
    CV_data = remove_delay(CV_data)
    CV_and_MS = synchronize([CV_data, MS_data_0])
    CV_and_MS['RE_vs_RHE'] = 0.553
    CV_and_MS['A_el'] = 0.2
    
    
    colors = {'M2':'b','M44':'r','M32':'k'}
    (ax1,ax2,ax3,) = plot_masses_and_I(CV_and_MS, colors=colors, tspan=CV_and_MS['tspan_2'], leg=0)
    (ax4,ax5,CV_and_MS_1) = plot_vs_potential(CV_and_MS, colors=colors, leg=0)
    
    
    
    
    