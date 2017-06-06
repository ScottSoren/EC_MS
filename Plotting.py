# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 19:07:45 2016
Most recently edited: 17B21

@author: scott
"""

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
import os

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from EC import sync_metadata
    from Data_Importing import import_folder
    from Combining import synchronize
    from Quantification import get_flux
    from Object_Files import lines_to_dictionary
    from Molecules import Molecule

    
else:                           #then we use relative import
    from .EC import sync_metadata
    from .Data_Importing import import_folder
    from .Combining import synchronize
    from .Quantification import get_flux
    from .Object_Files import lines_to_dictionary
    from .Molecules import Molecule
    
preferencedir = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'preferences' 
with open(preferencedir + os.sep + 'standard_colors.txt','r') as f:
    lines = f.readlines()
    standard_colors = lines_to_dictionary(lines)['standard colors']


def plot_vs_potential(CV_and_MS_0, 
                      colors={'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'},
                      tspan=0, RE_vs_RHE=None, A_el=None, 
                      ax1='new', ax2='new', ax=None, spec='k-',
                      overlay=0, logplot = [1,0], leg=1,
                      verbose=True, removebackground = None,
                      masses=None, mols=None, unit='nmol/s',
                      fig=None):
    '''
    This will plot current and select MS signals vs E_we, as is the 
    convention for cyclic voltammagrams. added 16I29
    '''
    if verbose:
        print('\n\nfunction \'plot_vs_potential\' at your service!\n')
    if type(logplot) is not list:
        logplot = [logplot, False]
    if removebackground is None:
        removebackground = not logplot[0]
    #prepare axes. This is ridiculous, by the way.
    CV_and_MS = CV_and_MS_0.copy() #17C01

    if ax == 'new':
        ax1 = 'new'
        ax2 = 'new'
    if ax1 != 'new':
        figure1 = ax1.figure
    elif ax2 != 'new':
        figure1 = ax2.figure
    else:
        if fig is None:
            figure1 = plt.figure()
        else:
            figure1 = fig
    if overlay:
        if ax1 == 'new':
            ax1 = figure1.add_subplot(111)
        if ax2 == 'new':
            ax2 = ax1.twinx()
    else:
        if ax1=='new':
            gs = gridspec.GridSpec(3, 1)
            gs.update(hspace=0.025)
            ax1 = plt.subplot(gs[0:2, 0])
            ax2 = plt.subplot(gs[2, 0])
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
    
    if ax2 is not None:
        #plot EC-lab data
        ax2.plot(V[I_plot],J[I_plot], spec)      
            #maybe I should use EC.plot_cycles to have different cycles be different colors. Or rewrite that code here.
        ax2.set_xlabel(V_str)
        ax2.set_ylabel(J_str)
    
    if ax1 is not None: #option of skipping an axis added 17C01
        #check if we're going to plot signals or fluxes:
        quantified = False      #added 17A07
        if mols is not None:
            colors = mols
            quantified = True
        elif ((type(colors) is dict and list(colors.keys())[0][0] == 'M') or
              (type(colors) is list and colors[0][0] == 'M' )):
            if masses is None:
                masses = colors
            else:
                colors = masses
        else:
            quantified = True
            mols = colors
        
        if type(colors) is list:
            c = colors.copy()
            colors = {}
            for m in c:
                if quantified:
                    mol = Molecule(m, verbose=False)
                    color = standard_colors[mol.primary]
                    colors[mol] = color
                else:
                    color = standard_colors[m]
                    colors[m] = color                
    
    #then do it.

        for (key, color) in colors.items():
            if quantified:
                (x,y) = get_flux(CV_and_MS, mol=key, tspan=tspan, removebackground=removebackground, 
                unit=unit, verbose=True)
                if type(key) is not str:
                    key = str(key) # in case key had been a Molecule object
                Y_str = key + '_' + unit
            else:   
                x_str = key + '-x'
                y_str = key + '-y'
                Y_str = key + '-Y'     #-Y will be a QMS signal interpreted to the EC-lab time variable.
                x = CV_and_MS[x_str]
                y = CV_and_MS[y_str]
            try:
                Y = np.interp(t, x, y)  #obs! np.interp has a has a different argument order than Matlab's interp1
            except ValueError:
                print('x ' + str(x) + '\ny ' + str(y) + '\nt ' + str(t))
            CV_and_MS[Y_str] = Y    #add the interpolated value to the dictionary for future use 
                                        #17C01: but not outside of this function.
            ax1.plot(V[I_plot], Y[I_plot], color, label=Y_str)
        if quantified:
            M_str = 'flux / [' + unit + ']'
        else:
            M_str = 'signal / [A]'
        #ax1.set_xlabel(V_str)
        ax1.set_xticks([])
        ax1.set_ylabel(M_str)
        if leg:
            ax1.legend()
    
    if verbose:
        print('\nfunction \'plot_vs_potential\' finished!\n\n')
        
        #parameter order of np.interp is different than Matlab's interp1
    return [ax1, ax2]    


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
        cols_1 = [[data_cols[i], data_cols[j]] for i,j in zip(I_axis_1[::2], I_axis_1[1::2])]        
            
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
        cols_2 = [[data_cols[i], data_cols[j]] for i,j in zip(I_axis_2[::2], I_axis_2[1::2])]

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
    '''something I used all the time back in the (Matlab) days.
        not sure I'll ever actually use it again though'''
    print(prompt + '\n... enter the indeces you\'re interested in, in order,' +
    'seperated by spaces, for example:\n>>>1 4 3')
    for nc, option in enumerate(options):
        print(str(nc) + '\t\t ' + options[nc])
    choice_string = input('\n')
    choices = choice_string.split(' ')
    choices = [int(choice) for choice in choices]
    return choices
    

def smooth_data(data_0, points=3, cols=None):
    '''
    Does a moving-average smoothing of data. I don't like it, but
    '''
    data = data_0.copy()
    if cols is None:
        cols = data['data_cols']
    for col in cols:
        x = data[col]
        c = np.array([1] * points) / points
        data[col] = np.convolve(x, c, mode='same')
    return data
        

def plot_signal(MS_data,
                masses = {'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'},
                tspan=0, ax='new', 
                logplot=True, saveit=False, leg=False, verbose=True):
    '''
    plots selected masses for a selected time range from MS data or EC_MS data
    Could probably be simplified a lot, to be the same length as plot_fluxes
    '''
    if verbose:
        print('\n\nfunction \'plot_masses\' at your service! \n Plotting from: ' + 
              MS_data['title'])

    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)    
    lines = {}
    if tspan == 0:                  #then use the range of overlap
        tspan = MS_data['tspan_2']  
    elif type(tspan) is str:
        tspan = MS_data[tspan]  
    if type(masses) is list:
        c = masses
        masses = {}
        for m in c:
            color = standard_colors[m]
            masses[m] = color

    for mass, color in masses.items():
        if verbose:
            print('plotting: ' + mass)
        x = MS_data[mass+'-x']
        y = MS_data[mass+'-y']
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
        lines[mass] = ax.plot(x, y, color, label = mass) 
        #as it is, lines is not actually used for anything         
    if leg:
        if type(leg) is not str:
            leg = 'lower right'
        ax1.legend(loc=leg)
    ax.set_xlabel('time / [s]')
    ax.set_ylabel('signal / [A]')           
    if logplot: 
        ax.set_yscale('log') 
    if verbose:
        print('function \'plot_masses\' finsihed! \n\n')
    return ax

def plot_masses(*args, **kwargs):
    print('plot_masses renamed plot_signal. Remember that next time!')
    return plot_signal(*args, **kwargs)
    
def plot_flux(MS_data, mols={'H2':'b', 'CH4':'r', 'C2H4':'g', 'O2':'k'},
            tspan='tspan_2', ax='new', removebackground=True, A_el=None,
            logplot=True, leg=False, unit='nmol/s', verbose=True):
    '''
    Plots the molecular flux to QMS in nmol/s for each of the molecules in
    'fluxes.keys()', using the primary mass and the F_cal value read from
    the molecule's text file, with the plot specs from 'fluxes.values()'
    '''
    if verbose:
        print('\n\nfunction \'plot_flux\' at your service!\n')
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)  
    if type(tspan) is str:
        tspan = MS_data[tspan]
    if type(mols) is list:
        c = mols
        mols = {}
        for m in c:
            mol = Molecule(m, verbose=False)
            color = standard_colors[mol.primary]
            mols[mol] = color
        
    for (mol, color) in mols.items():
        [x,y] = get_flux(MS_data, mol, unit=unit, verbose=verbose, tspan=tspan)
        '''  17A28: this is now taken care of in the line above.
        if tspan is not None:
            if verbose:
                print('cutting ' + str(mol) + ' at ' + str(tspan))
            x,y = cut(x, y, tspan)
        '''
        if removebackground:
            try:
                y = y - 0.99 * min(y) #0.99 to avoid issues when log plotting.
            except ValueError:
                print(y)
        if type(mol) is str:
            l = mol
        else:
            l = mol.name
        ax.plot(x, y, color, label=l)
    if leg:
        if type(leg) is not str:
            leg = 'lower right'
        ax.legend(loc=leg)
    ax.set_xlabel('time / [s]')
    ylabel = 'flux / [' + unit + ']'

    ax.set_ylabel(ylabel)
    if logplot:
        ax.set_yscale('log')
    
    if verbose:
        print('\nfunction \'plot_flux\' finished!\n\n')    
    return ax    

    
def plot_experiment(EC_and_MS,
                    colors={'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'},
                    tspan=None, overlay=False, logplot=[True,False], verbose=True,   
                    plotpotential=True, plotcurrent=True, ax='new',
                    RE_vs_RHE=None, A_el=None, removebackground=True,
                    saveit=False, title=None, leg=False, unit='nmol/s',
                    masses=None, mols=None, #mols will overide masses will overide colors
                    potentialcolor='k', currentcolor='r', 
                    potentiallabel=None, currentlabel=None,
                    fig=None,
                    ): 
    '''
    this plots signals or fluxes on one axis and current and potential on other axesaxis
    '''
    
    if verbose:
        print('\n\nfunction \'plot_experiment\' at your service!\n Plotting from: ' + EC_and_MS['title'])
    if ax == 'new':
        if fig is None:
            figure1 = plt.figure()
        else:
            figure1 = fig
            plt.figure(figure1.number)
            print('plot_expeiriment using ' + str(fig))
        if overlay:
            ax = [figure1.add_subplot(111)]
            ax += [ax[0].twinx()]                     
        else:
            gs = gridspec.GridSpec(3, 1)
            gs.update(hspace=0.025)
            ax = [plt.subplot(gs[0:2, 0])]
            ax += [plt.subplot(gs[2, 0])]
            if plotcurrent and plotpotential:
                ax += [ax[1].twinx()]
        
    if tspan is None:                  #then use the whole range of overlap
        tspan = EC_and_MS['tspan_2']
    if type(logplot) is not list:
        logplot = [logplot, False]

    V_str, J_str = sync_metadata(EC_and_MS, RE_vs_RHE=RE_vs_RHE, A_el=A_el) #added 16J27
    A_el = EC_and_MS['A_el']

    quantified = False      #added 16L15
    if mols is not None:
        quantified = True
    elif ((type(colors) is dict and list(colors.keys())[0][0] == 'M') or
          (type(colors) is list and colors[0][0] == 'M' )):
        if masses is None:
            masses = colors
    else:
        quantified = True
        mols = colors
    
    if quantified:
        plot_flux(EC_and_MS, mols=mols, tspan=tspan, A_el=A_el,
                  ax=ax[0], leg=leg, logplot=logplot[0], unit=unit, 
                  removebackground=removebackground, verbose=verbose)
    else:
        plot_signal(EC_and_MS, masses=masses, tspan=tspan,
                    ax=ax[0], leg=leg, logplot=logplot[0], verbose=verbose)
    if not overlay:
        ax[0].set_xlabel('')
        ax[0].xaxis.tick_top()  
    
    if title is not None:
            plt.title(title)
    
    t = EC_and_MS['time/s']
        
    V = EC_and_MS[V_str]
    J = EC_and_MS[J_str]      
    
    if tspan is not None:
        I_keep = [I for (I, t_I) in enumerate(t) if tspan[0]<t_I and t_I<tspan[1]]
        t = t[I_keep]
        V = V[I_keep]
        J = J[I_keep]

    i_ax = 1
    if plotpotential:
        ax[i_ax].plot(t, V, color=potentialcolor, label=potentiallabel)
        ax[i_ax].set_ylabel(V_str)
        if len(logplot) >2:
            if logplot[2]:
                ax[i_ax].set_yscale('log')
        xlim = ax[i_ax-1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        i_ax += 1
        
    if plotcurrent:
        ax[i_ax].plot(t, J, currentcolor, label=currentlabel)
        ax[i_ax].set_ylabel(J_str)
        ax[i_ax].set_xlabel('time / [s]')
        xlim = ax[i_ax-1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        if logplot[1]: 
            ax[i_ax].set_yscale('log')
    if plotcurrent or plotpotential:
        ax[1].set_xlabel('time / [s]')
    
    if saveit:
        if title == 'default':
            title == EC_and_MS['title'] + '.png'
        figure1.savefig(title)
        
    if verbose:
        print('function \'plot_experiment\' finished!\n\n')
    
    return ax
    
def plot_masses_and_I(*args, **kwargs):
    print('plot_masses_and_I renamed plot_experiment. Remember that next time!')
    return plot_experiment(*args, **kwargs)

def plot_folder(folder_name, 
                colors={'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'}, 
                RE_vs_RHE=None, A_el=None):
    '''
    Plots an EC and MS data from an entire folder, generally corresponding to
    a full day of measurements on one sample.
    Will probably only be used to get an overview.
    Could add text showing starts of the data files
    '''
    Datasets = import_folder(folder_name)
    Combined_data = synchronize(Datasets, t_zero='first')
    sync_metadata(Combined_data, RE_vs_RHE, A_el)
    return plot_experiment(Combined_data, colors=colors)


def plot_datapoints(integrals, colors, ax='new', label='', X=None, X_str='V',
                    logplot=True, specs={}, Xrange=None):
    '''
    integrals will most often come from functino 'get_datapoitns' in module
    Integrate_Signals
    '''
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    if X is None:
        X = integrals[X_str]
        
    for (quantity, color) in colors.items(): 
        # Here I just assme they've organized the stuff right to start with.
        # I could alternately use the more intricate checks demonstrated in
        # DataPoints.plot_errorbars_y
        value = integrals[quantity]
        if type(Xrange) is dict:
            Xrange_val = Xrange[quantity]
        else:
            Xrange_val = Xrange 
        if type(color) is dict:
            plot_datapoints(value, color, ax=ax, logplot=logplot,
                            label=label+quantity+'_', X=X, Xrange=Xrange_val, specs=specs)
        else:
            if type(color) is tuple: #note a list can be a color in rbg
                spec = color[0]
                color = color[1]
                if 'markersize' not in specs:
                    specs['markersize'] = 5
            else:
                spec = '.'
                if 'markersize' not in specs:
                    specs['markersize'] = 15
            #print(quantity + '\n\tvalue=' + str(value) + 
            #        '\n\tcolor=' + str(color) + '\n\tV=' + str(V))
            #print(quantity + ' ' + str(color))
            if Xrange is not None:
                I_keep = np.array([I for (I, X_I) in enumerate(X) if 
                          Xrange_val[0] <= float(np.round(X_I,2)) <= Xrange_val[1]])
                X_plot = np.array(X)[I_keep]
                value_plot = np.array(value)[I_keep]
                #there was a mindnumbing case of linking here.
                #tried fix it with .copy(), but new variable names needed.
            else:
                X_plot = X
                value_plot = value
            ax.plot(X_plot, value_plot, spec, 
                    color=color, label=label+quantity, **specs, )
    if logplot:
        ax.set_yscale('log')
    return ax

    
if __name__ == '__main__':
    import os
    from Data_Importing import import_data
    from EC import select_cycles, remove_delay
    
    plt.close('all')
    
    importrawdata = 1
    if importrawdata:
        default_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) 
        CV_data_0 = import_data(default_directory + os.sep, # + '18_CO_dose_and_strip_C01.mpt', data_type='EC')
                                data_type='EC')        
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
    
    
    
    
    