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
#from mpl_toolkits.axes_grid1 import make_axes_locatable

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from EC import sync_metadata, select_cycles
    from Data_Importing import import_folder
    from Combining import synchronize
    from Quantification import get_flux, get_signal
    from Object_Files import lines_to_dictionary
    from Molecules import Molecule

    
else:                           #then we use relative import
    from .EC import sync_metadata, select_cycles
    from .Data_Importing import import_folder
    from .Combining import synchronize
    from .Quantification import get_flux, get_signal
    from .Object_Files import lines_to_dictionary
    from .Molecules import Molecule
    
preferencedir = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'preferences' 
with open(preferencedir + os.sep + 'standard_colors.txt','r') as f:
    lines = f.readlines()
    standard_colors = lines_to_dictionary(lines)['standard colors']

def get_standard_colors():
    return standard_colors


def plot_vs_potential(CV_and_MS_0, 
                      colors={'M2':'b','M4':'m','M18':'y','M28':'0.5','M32':'k'},
                      tspan=0, RE_vs_RHE=None, A_el=None, cycles='all',
                      ax1='new', ax2='new', ax=None, #spec='k-',
                      overlay=0, logplot = [1,0], leg=False,
                      verbose=True, removebackground = None,
                      masses=None, mols=None, unit=None,
                      fig=None, spec={}):
    '''
    This will plot current and select MS signals vs E_we, as is the 
    convention for cyclic voltammagrams. added 16I29
    
    #there's a lot of code here that's identical to plot_experiment. Consider
    #having another function for e.g. processing these inputs.
    '''
    if verbose:
        print('\n\nfunction \'plot_vs_potential\' at your service!\n')
    if type(logplot) is not list:
        logplot = [logplot, False]
    if removebackground is None:
        removebackground = not logplot[0]
    #prepare axes. This is ridiculous, by the way.
    CV_and_MS = CV_and_MS_0.copy() #17C01
    if not cycles == 'all':
        CV_and_MS = select_cycles(CV_and_MS, cycles, verbose=verbose)

    if ax == 'new':
        ax1 = 'new'
        ax2 = 'new'
    elif ax is not None:
        ax1 = ax[0]
        ax2 = ax[1]
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
            #gs.update(hspace=0.025)
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
        tspan = CV_and_MS['tspan']
    I_plot = np.array([i for (i,t_i) in enumerate(t) if tspan[0]<t_i and t_i<tspan[1]])
    
    if ax2 is not None:
        #plot EC-lab data
        ec_spec = spec.copy()
        if 'color' not in ec_spec.keys():
            ec_spec['color'] = 'k'
        ax2.plot(V[I_plot],J[I_plot], **ec_spec)      
            #maybe I should use EC.plot_cycles to have different cycles be different colors. Or rewrite that code here.
        ax2.set_xlabel(V_str)
        ax2.set_ylabel(J_str)
    
#    print(masses)
    if ax1 is not None: #option of skipping an axis added 17C01
        #check if we're going to plot signals or fluxes:
        quantified = False      #added 16L15
        if mols is not None:
            quantified = True
            colors = mols  #added 17H11  
        elif masses is not None:
#            print('masses specified')
            quantified = False
            colors = masses
        elif ((type(colors) is dict and list(colors.keys())[0][0] == 'M') or
              (type(colors) is list and type(colors[0]) is str and colors[0][0] == 'M' ) or
              (type(colors) is str and colors[0]=='M')):
            print('uncalibrated data to be plotted.')
            masses = colors    
            colors = masses
        else:
            quantified = True
            mols = colors
        if type(colors) is not list and type(colors) is not dict:
            colors = [colors]        
        
#        print(type(colors))
        
        if unit is None:
            if quantified:
                unit = 'pmol/s'
            else:
                unit = 'nA'
        if type(colors) is list:
            c = colors.copy()
            colors = {}
            for m in c:
                print(str(m))
                if quantified:
                    if type(m) is str:
                        mol = Molecule(m, verbose=False)
                    else:
                        mol = m  
                    color = standard_colors[mol.primary]
                    colors[mol] = color
                else:
                    color = standard_colors[m]
                    colors[m] = color                
    
    #then do it.

        for (key, color) in colors.items():
            if quantified:
                (x,y) = get_flux(CV_and_MS, mol=key, tspan=tspan, 
                removebackground=removebackground, unit=unit, verbose=True)
                if type(key) is not str:
                    key = str(key) # in case key had been a Molecule object
                Y_str = key + '_' + unit
            else:   

                Y_str = key + '_' + unit    #        
                x, y = get_signal(CV_and_MS, mass=key, tspan=tspan, unit=unit)

            try:
                Y = np.interp(t, x, y)  #obs! np.interp has a has a different argument order than Matlab's interp1
            except ValueError:
                print('x ' + str(x) + '\ny ' + str(y) + '\nt ' + str(t))
            CV_and_MS[Y_str] = Y    #add the interpolated value to the dictionary for future use 
                                        #17C01: but not outside of this function.
            ms_spec = spec.copy()
            if 'color' not in ms_spec.keys():
                ms_spec['color'] = color
            ax1.plot(V[I_plot], Y[I_plot], label=Y_str, **ms_spec)
        if quantified:
            M_str = 'cal. signal / [' + unit + ']'
        else:
            M_str = 'MS signal / [' + unit + ']'
        #ax1.set_xlabel(V_str)
        ax1.xaxis.tick_top()
        ax1.set_ylabel(M_str)
        if leg:
            ax1.legend()
    
    if verbose:
        print('\nfunction \'plot_vs_potential\' finished!\n\n')
    
    for ax in [ax1, ax2]:    
        ax.tick_params(axis='both', direction='in') #17K28  
        
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
        
    axes_1.set_xlabel('time / [s]')
    axes_1.set_ylabel('log(MS signal / [a.u.])')
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
    axes_2.set_ylabel('current / [mA]')
    axes_2.set_xlabel('time / [s]')
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
    

def smooth_data(data_0, points=3, cols=None, verbose=True):
    '''
    Does a moving-average smoothing of data. I don't like it, but
    experencing problems 17G26
    '''
    data = data_0.copy()
    if cols is None:
        cols = data['data_cols']
    for col in cols:
        if verbose:
            print('smoothening \'' + col + '\' with a ' + str(points) + '-point moving average')
        x = data[col]
        c = np.array([1] * points) / points
        #print(str(len(c)))
        #data[col] = 0 #somehow this doesn't come through, 17G25
        data[col] = np.convolve(x, c, mode='same')
        #print('len = ' + str(len(x)))
        #x = None

    data['test'] = None #This does make it through. 
    #data['U vs RHE / [V]'] = None #this doesn't make it through either.
    return data    

def plot_signal(MS_data,
                masses = {'M2':'b','M4':'r','M18':'0.5','M28':'g','M32':'k'},
                tspan=None, ax='new', unit='nA',
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
    #note, tspan is processed in get_signal, and not here!
    if type(masses) is str:
        masses = [masses]
    if type(masses) is list:
        c = masses
        masses = {}
        for m in c:
            color = standard_colors[m]
            masses[m] = color

    for mass, color in masses.items():
        if verbose:
            print('plotting: ' + mass)
        try:
            x, y = get_signal(MS_data, mass, unit=unit, verbose=verbose, tspan=tspan)
        except KeyError:
            print('Can\'t get signal for ' + str(mass))
            continue
        lines[mass] = ax.plot(x, y, color, label = mass) 
        #as it is, lines is not actually used for anything         
    if leg:
        if type(leg) is not str:
            leg = 'lower right'
        ax1.legend(loc=leg)
    ax.set_xlabel('time / [s]')
    ax.set_ylabel('MS signal / [' + unit + ']')           
    if logplot: 
        ax.set_yscale('log') 
    ax.tick_params(axis='both', direction='in') #17K28  
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
    #note, tspan is processed in get_flux, and not here!
    
    if type(mols) is not list and type(mols) is not dict:
        #then it's probably just one molecule object.
        mols = [mols]
    
    if type(mols) is list:
        c = mols
        mols = {}
        for m in c:
            if type(m) is str:
                mol = Molecule(m, verbose=False)
            else:
                mol = m  #this function should accept a list of Molecule instances!
            color = standard_colors[mol.primary]
            mols[mol] = color
    
    for (mol, color) in mols.items():
        try:
            [x,y] = get_flux(MS_data, mol, unit=unit, verbose=verbose, tspan=tspan)
        except KeyError:
            print('Can\'t get signal for ' + str(mol))
            continue
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
    ylabel = 'cal. signal / [' + unit + ']'

    ax.set_ylabel(ylabel)
    if logplot:
        ax.set_yscale('log')
        
    ax.tick_params(axis='both', direction='in') #17K28  
    
    if verbose:
        print('\nfunction \'plot_flux\' finished!\n\n')    
    return ax    

    
def plot_experiment(EC_and_MS,
                    colors={'M2':'b','M4':'m','M18':'y','M28':'0.5','M32':'k'},
                    tspan=None, overlay=False, logplot=[True,False], verbose=True,   
                    plotpotential=True, plotcurrent=True, ax='new',
                    RE_vs_RHE=None, A_el=None, removebackground=True,
                    saveit=False, title=None, leg=False, unit='pmol/s',
                    masses=None, mols=None, #mols will overide masses will overide colors
                    V_color='k', J_color='r', V_label=None, J_label=None,
                    fig=None, t_str=None, J_str=None, V_str=None
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
            #gs.update(hspace=0.025)
            #gs.update(hspace=0.05)
            ax = [plt.subplot(gs[0:2, 0])]
            ax += [plt.subplot(gs[2, 0])]
            if plotcurrent and plotpotential:
                ax += [ax[1].twinx()]
        
    if tspan is None:                  #then use the whole range of overlap
        if 'tspan_EC' in EC_and_MS:
            tspan = EC_and_MS['tspan_EC']
        else:                              
            tspan = EC_and_MS['tspan'] #changed from 'tspan_2' 17H09
    if type(tspan) is str and not tspan=='all':
        tspan = EC_and_MS[tspan]
    if type(logplot) is not list:
        logplot = [logplot, False]
    
    if t_str is None:
        t_str = 'time/s'
    if V_str is None or J_str is None or RE_vs_RHE is not None or A_el is not None: 
        V_str_0, J_str_0 = sync_metadata(EC_and_MS, RE_vs_RHE=RE_vs_RHE, A_el=A_el, verbose=verbose) 
        #added 16J27... problem caught 17G26, fixed in sync_metadata
    if V_str is None: #this way I can get it to plot something other than V and J.
        V_str = V_str_0
    if J_str is None:
        J_str = J_str_0
    
    A_el = EC_and_MS['A_el']

    quantified = False      #added 16L15
    #print(type(colors))
    #if type(colors) is list and type(colors[0]) is not str:
    #    print(type(colors[0]))
    if mols is not None:
        quantified = True
    elif ((type(colors) is dict and list(colors.keys())[0][0] == 'M') or
          (type(colors) is list and type(colors[0]) is str and colors[0][0] == 'M' ) or
          (type(colors) is str and colors[0]=='M')):
        print('uncalibrated data to be plotted.')
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
    
    ax[0].set_xlim(tspan)
    
    if title is not None:
            plt.title(title)
    
    try:
        t = EC_and_MS[t_str]       
    except KeyError:
        print('data doesn\'t contain \'' + str(t_str) + '\', i.e. t_str. Can\'t plot EC data.')
        plotpotential = False
        plotcurrent = False
    try:
        V = EC_and_MS[V_str]
    except KeyError:
        print('data doesn\'t contain \'' + str(V_str) + '\', i.e. V_str. Can\'t plot that data.')
        plotpotential = False      
    try:
        J = EC_and_MS[J_str]      
    except KeyError:
        print('data doesn\'t contain \'' + str(J_str) + '\', i.e. J_str. Can\'t plot that data.')
        plotcurrent = False     
        
        # to check if I have problems in my dataset
#    print('len(t) = ' + str(len(t)) + 
#          '\nlen(V) = ' + str(len(V)) + 
#          '\nlen(J) = ' + str(len(J)))
    
    if tspan is not 'all' and plotcurrent or plotpotential:
        I_keep = [I for (I, t_I) in enumerate(t) if tspan[0]<t_I and t_I<tspan[1]]
        t = t[I_keep]
        if plotpotential:
            V = V[I_keep]
        if plotcurrent:
            J = J[I_keep]

    i_ax = 1
    if plotpotential:
        ax[i_ax].plot(t, V, color=V_color, label=V_label)
        ax[i_ax].set_ylabel(V_str)
        if len(logplot) >2:
            if logplot[2]:
                ax[i_ax].set_yscale('log')
        xlim = ax[i_ax-1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        ax[i_ax].yaxis.label.set_color(V_color)
        ax[i_ax].tick_params(axis='y', colors=V_color)
        ax[i_ax].spines['left'].set_color(V_color)
        ax[i_ax].tick_params(axis='both', direction='in') #17K28  
        i_ax += 1
        
    if plotcurrent:
        ax[i_ax].plot(t, J, color=J_color, label=J_label)
        ax[i_ax].set_ylabel(J_str)
        ax[i_ax].set_xlabel('time / [s]')
        xlim = ax[i_ax-1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        if logplot[1]: 
            ax[i_ax].set_yscale('log')
        ax[i_ax].yaxis.label.set_color(J_color)
        ax[i_ax].tick_params(axis='y', colors=J_color)
        if i_ax == 2:
            ax[i_ax].spines['right'].set_color(J_color)
        else:
            ax[i_ax].spines['left'].set_color(J_color)
        ax[i_ax].tick_params(axis='both', direction='in') #17K28  
        
    if plotcurrent or plotpotential:
        ax[1].set_xlabel('time / [s]')
        ax[1].set_xlim(tspan)
    
    if saveit:
        if title == 'default':
            title == EC_and_MS['title'] + '.png'
        figure1.savefig(title)
        
    if verbose:
        print('function \'plot_experiment\' finished!\n\n')
    
    return ax
    
def plot_masses_and_I(*args, **kwargs):
    print('\n\n\'plot_masses_and_I\' has been renamed \'plot_experiment\'. Remember that next time!')
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




def plot_operation(cc=None, t=None, j=None, z=None, tspan=None, results=None,
                   plot_type='heat', ax='new', colormap='plasma', aspect='auto', 
                   unit='pmol/s', dimensions=None, verbose=True):
    if verbose:
        print('\n\nfunction \'plot_operation\' at your service!\n')
    # and plot! 
    if type(cc) is dict and results is None:
        results = cc
        cc = None
    if results is None:
        results = {} #just so I don't get an error later
    if cc is None:
        cc = results['cc']
    if t is None:
        if 't' in results:
            t = results['t']
        elif 'x' in results:
            t = results['x']
        else:
            t = np.linspace(0, 1, np.size(cc, axis=0))
    if z is None:
        if 'z' in results:
            z = results['z']
        elif 'y' in results:
            z = results['y']
        else:
            z = np.linspace(0, 1, np.size(cc, axis=1))
    if j is None:
        if j in results:
            j = results['j']
        else:
            j = cc[0, :]   
    if dimensions is None:
        if 'dimensions' in results:
            dimensions = results['dimensions']
        else:
            dimensions = 'tz'
    
    
    if tspan is None:
        tspan = [t[0], t[-1]]
    if plot_type == 'flux':   
        if ax == 'new':
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
        else:
            ax1 = ax #making a heat map will only work with a new axis.
        ax1.plot(t, j, label='simulated flux')
        ax1.set_xlabel('time / [s]')
        ax1.set_ylabel('flux / [' + unit + ']')
        axes = ax1
        
    elif plot_type == 'heat' or plot_type == 'both':
        if ax == 'new':
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
        elif type(ax) is list: 
            ax1 = ax[0]
        else:
            ax1 = ax
        
        #t_mesh, x_mesh = np.meshgrid(t,x)
        #img = ax1.contourf(t_mesh, x_mesh*1e6, np.transpose(cc,[1,0]), cmap='Spectral', 
        #                   levels=np.linspace(np.min(cc),np.max(cc),100))
        
        # imshow objects seem more versatile than contourf for some reason.
        
        trange = [min(t), max(t)]
        if dimensions[0] == 'x':
            trange = [t*1e3 for t in trange] # m to mm
        zrange = [min(z*1e6), max(z*1e6)]
        
        img = ax1.imshow(np.transpose(cc,[1,0]), 
                         extent=trange[:] + zrange[:],  #have to be lists here!
                         aspect=aspect, origin='lower',
                         cmap = colormap)        
        
#        divider = make_axes_locatable(ax1)
#        cax = divider.append_axes("right", size="5%", pad=0.05)
#https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph

        cbar = plt.colorbar(img, ax=ax1)
        cbar.set_label('concentration / [mM]')
        if dimensions[0] == 't':
            ax1.set_xlabel('time / [s]')
        elif dimensions[0] == 'x':
            ax1.set_xlabel('position / [mm]')
        ax1.set_ylabel('position / [um]')
        
#        print('plot_type = ' + plot_type)
        if plot_type == 'both':
            if type(ax) is list:
                ax2 = ax[1]
            else:
                ax2 = ax1.twinx()
            ax2.set_ylabel('flux / [' + unit + ']')
            ax2.plot(t, j, 'k-')
            cbar.remove()
            ax3 = img.figure.add_axes([0.85, 0.1, 0.03, 0.8])
            cbar = plt.colorbar(img, cax=ax3)
            cbar.set_label('concentration / [mM]')
            ax1.set_xlim(tspan)
            print('returning three axes!')
            axes = [ax1, ax2, ax3]
        else:
            axes = [ax1, cbar]
            
    if verbose:
        print('\nfunction \'plot_operation\' finished!\n\n')
    return axes

def set_figparams(figwidth=8,aspect=4/3,fontsize=7,figpad=0.15,figgap=0.08):
    import matplotlib as mpl
    
    #figwidth=8  #figwidth in cm width 20.32cm = 8inches being standard and thesis textwidth being 12.
    #aspect=4/3  #standard is 4/3

    #fontsize=7  #standard is 12.0pt, thesis is 10.0pt and footnotesize is 8.0pt and lower case seems to be 2/3 of full font size, which makes 7pt "nice" for thesis plotting
    
    realfigwidth=20*(fontsize/12)*1.2  #a factor 1.2 makes all lines "equaly thick" - make this 1 for figwidth=4cm (sniff2fig6)
    #figsize=[20.32,15.24]
    figsize=[realfigwidth,realfigwidth/aspect]
    
    mpl.rc('font', size=fontsize*(realfigwidth/figwidth))
    
    mpl.rc('mathtext', fontset='custom',
                       rm='Helvetica',
                       #it='Helvetica:italic',
                       #bf='Helvetica:bold',
                       )
    
    mpl.rc('figure', figsize=[figsize[0]/2.54,figsize[1]/2.54],
                     dpi=100*2.54*figwidth/realfigwidth
                     )
    
    #figpad=0.14  #fraction of figure size
    #figgap=0.08  #fraction of figure size
    mpl.rc('figure.subplot', left=figpad,
                             right=1-figpad,
                             bottom=figpad,
                             top=1-figpad,
                             hspace=figgap,
                             )
    mpl.rc('xtick', labelsize='small') 
    mpl.rc('ytick', labelsize='small')
    
    #mpl.rc('axes', labelweight='medium')
    
    mpl.rc('savefig', dpi=250*2.54*figwidth/realfigwidth)

    
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
    
    
    
    
    
