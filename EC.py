

# -*- coding: utf-8 -*-
"""
created 16I15
last edited 16J27

@author: Scott

functions for use on EC-lab data
"""

#make python2-compatible
from __future__ import print_function
from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
import os

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from Combining import timestamp_to_seconds, is_time, cut
else:                           #then we use relative import
    from .Combining import timestamp_to_seconds, is_time, cut


def select_cycles(EC_data_0, cycles=1, t_zero=None, verbose=True):
    ''' 
    This function selects one or more cycles from EC_data_0.
    Use this before synchronizing!
    Works for both CA and CV
    #changed 16L22 to work on EC_and_MS data
    '''
    if verbose:
        print('\nSelecting cycles ' + str(cycles) + ' from \'' + EC_data_0['title'] + '\'\n')
    
    EC_data = EC_data_0.copy()
    
    #it looks like I actually want Ns for CA's
    if 'Ns' in EC_data['data_cols']:
        cycle_numbers = EC_data['Ns']
    else:
        cycle_numbers = EC_data['cycle number']
    
    N = len(cycle_numbers)
    if type(cycles)==int:
        cycles = [cycles]
    I_keep = np.array([I for I in range(N) if cycle_numbers[I] in cycles])
    #list comprehension is awesome.

    for col in EC_data['data_cols']:
        try:
            if not (col[0] == 'M' and col[-2:] in ['-x', '-y']):  
                #then we're dealing with EC data
                try:
                    EC_data[col] = EC_data[col].copy()[I_keep]
                except KeyError:
                    print('hm... \'' + col + '\' in EC_data[data_cols] but not in EC_data')
        except IndexError:
            print('trouble selecting cycle ' + str(cycles) + ' of ' + col + '\n' +
                    'type(I_keep) = ' + str(type(I_keep)))
    t0 = timestamp_to_seconds(EC_data['timestamp'])
    tspan_2 = np.array([min(EC_data['time/s']), max(EC_data['time/s'])])
    EC_data['tspan_2'] = tspan_2
    EC_data['tspan'] = tspan_2 + t0
    EC_data['data_type'] += ' selected'   
    
    for col in EC_data['data_cols']:
        if col[1] == 'M' and col[-2:] == '-x': #then we've got a QMS time variable
            y_col = col[:-2] + '-y'
            EC_data[col], EC_data[y_col] = cut(EC_data[col], EC_data[y_col], tspan_2)       

    if t_zero is 'start':
        for col in EC_data['data_cols']:
            if is_time(col):
                EC_data[col] = EC_data[col] - tspan_2[0]
                EC_data['tspan_2'] = tspan_2 - tspan_2[0]
    
    return EC_data


def remove_delay(CV_data):
    '''
    Gets rid of the delay at the beginning of .mpt files before it actually starts 
    cycling. This is not seen cycle_number, but in control changes, which goes to 0 for 
    the first time right as the cycle starts... I think. 16I29
    '''
    control = CV_data['control changes']
    
    I_start = np.where(control == 0)[0][0]

    for col in CV_data['data_cols']:
        CV_data[col] = CV_data[col][I_start:-1]

    return CV_data
        

def CV_difference(cycles_data, redox=1, Vspan=[0.5, 1.0], 
                  ax=None, color='g', verbose=1):
    '''
    This will calculate the difference in area between two cycles in a CV, 
    written for CO stripping 16J26. If ax is given, the difference will be
    filled in with color.
    '''
    if verbose:
        print('\n\nfunction \'CV_difference\' at your service!\n')  
    
    if redox == 'ox':
        redox = [1]
    elif redox == 'red':
        redox = [0]
    elif type(redox) is int:
        redox = [redox]
    elif redox is None:
        redox = [0,1]

    Vs = []   
    Js = []
    Q = []
    JV = []
    ts = []
    for cycle_data in cycles_data:
        V_str, J_str = sync_metadata(cycle_data)
        V = cycle_data[V_str]
        J = cycle_data[J_str]
        t = cycle_data['time/s']
        
        ro = cycle_data['ox/red']
        q = cycle_data['(Q-Qo)/C']
        I_keep = [I for (I, (V, ro)) in enumerate(zip(V, ro)) if 
                    Vspan[0] < V < Vspan[1] and ro in redox]
        
        V = V[I_keep]        
        J = J[I_keep]
        t = t[I_keep]
        
        Vs += [V]
        Js += [J]
        ts += [t]
        
        Q += [q[I_keep[-1]] - q[I_keep[0]]]       
        JV += [np.trapz(J, V)]
        
    dQ = Q[0] - Q[1] 
    dJV = JV[0] - JV[1]
    
    if verbose:
        try:
            A_el = cycle_data['A_el']
        except KeyError:
            A_el = 1
            print('didn''t find A_el.')
        print('difference in charge passed: a = ' + str(dQ) + ' C\n' + 
                'difference in CV area: b = ' + str(dJV) + ' V*mA/cm^2\n' + 
                'scan rate: b/a*A_el = ' + str(dJV / dQ * A_el) + ' mV/s') 
    
    if len(Vs[0]) != len(Vs[1]):  #then we'll have to interpolate
        if 1 in redox:
            Js[1] = np.interp(Vs[0], Vs[1], Js[1])
            V_avg = Vs[0]
        else:    
            Js[1] = np.interp(-Vs[0], -Vs[1], Js[1])
    else:
        V_avg = (Vs[0] + Vs[1]) / 2
    J_diff = Js[0] - Js[1] #note this is all optimized for CO stripping
    
    if ax:
        ax.fill_between(V_avg, Js[0], Js[1], where=Js[0]>Js[1],
                        facecolor=color, interpolate=True)
    if verbose:
        print('\nfunction \'CV_difference\' finished!\n\n')
    
    return dQ, [ts, V_avg, J_diff]

    
def smooth_pulses(CA_Data_0, verbose=1):
    '''
    This function turns the CA data into a square wave by averaging the
    potential over the duration of a pulse (where it should be constant).
    Useful when noise makes the figures look ugly dispite otherwise good data.
    If you have to use this function, though, I would say the results are not 
    publication-ready.
    '''
    if verbose:
        print('\n\nfunction \'smooth_pulses\' at your service!')
    CA_Data = CA_Data_0.copy()
    cycle_numbers = CA_Data['Ns']
    cycles = np.unique(cycle_numbers)
    for c in cycles:
        I_cycle = np.array([i for (i,cycle) in enumerate(cycle_numbers) if cycle==c])
        V_avg = np.average(CA_Data['Ewe/V'][I_cycle])
        CA_Data['Ewe/V'][I_cycle] = V_avg
    if verbose:
        print('function \'smooth_pulses\' finished!\n\n')
    return CA_Data


def plot_vs_time(EC_data, axes='new', y_strings='default', colors=None,
                 RE_vs_RHE=None, A_el=None, verbose=1):
    
    if verbose:
        print('\n\nfunction \'plot_vs_time\' at your service!')    
    
    V_str, J_str = sync_metadata(EC_data, RE_vs_RHE, A_el)
    if y_strings == 'default':
        y_strings = V_str, J_str
    
    t_str = 'time/s'
    t = EC_data[t_str]
    
    if colors is None:
        colors = ['k'] * len(y_strings)
    
    if axes == 'new':
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax2 = ax1.twinx()
        axes = [ax1, ax2]
    
    for (ax, y_str, color) in zip(axes, y_strings, colors):
        try:
            y = EC_data[y_str]
        except KeyError:
            print('Can''t find ' + y_str + '. skipping that one.')
            continue
        
        ax.plot(t, y, color, label=y_str)
        ax.set_xlabel(t_str)        
        ax.set_ylabel(y_str)
    
    if verbose:
        print('function \'plot_vs_time\' finished!\n\n')  
        
    return axes


def sync_metadata(EC_data, RE_vs_RHE=None, A_el=None, verbose=1):
    '''
    Deal with all the annoying RE and J vs I vs <I> stuff once and for all here.
    After this has been called, all plotting methods need only to utilize
    EC_data['V_str'] and EC_data['J_str']
    '''    
    if verbose:
        print('\nsyncing metadata for ' + EC_data['title'] + '\n')
        
    if RE_vs_RHE is not None:
        EC_data['RE_vs_RHE'] = RE_vs_RHE
    elif 'RE_vs_RHE' in EC_data:
        RE_vs_RHE = EC_data['RE_vs_RHE']
    else:
        EC_data['RE_vs_RHE'] = None
    
    if RE_vs_RHE is None:
        V_str = 'Ewe/V'
    else:
        V_str = 'E vs RHE / [V]'
        EC_data[V_str] = EC_data['Ewe/V'] + RE_vs_RHE
    
    if A_el is not None:
        EC_data['A_el'] = A_el
    elif 'A_el' in EC_data:
        A_el = EC_data['A_el']
    else:
        EC_data['A_el'] = None
        
    I_str = [s for s in ['I/mA', '<I>/mA'] if s in EC_data['data_cols']][0]     
    if A_el is None:
        J_str = I_str
    else:
        J_str = 'J /[mA/cm^2]'
        EC_data[J_str] = EC_data[I_str] / A_el
 
    EC_data['V_str'] = V_str   
    EC_data['J_str'] = J_str
    EC_data['I_str'] = I_str
    
    for col in [V_str, J_str]:
        if col not in EC_data['data_cols']:
            EC_data['data_cols'] += [col]
            if verbose:
                print('added ' + col + ' to data_cols')
        
    return V_str, J_str
    

def plot_CV_cycles(CV_data, cycles=[0], RE_vs_RHE=None, A_el=None, ax='new',
                   saveit=0, title='default', leg=0, verbose=1, colors=None):
    '''
    plots a subset of cycles in a CV
    '''
    if verbose:
        print('\n\nfunction \'plot_CV_cycles\' at your service!')
    
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111) 
    
    if type(cycles)==int:
        cycles = [cycles]  
        
    V_str, J_str = sync_metadata(CV_data, RE_vs_RHE, A_el) #added 16J26
    
    data_to_return = [] 
    for n, cycle in enumerate(cycles):

        cycle_data = select_cycles(CV_data, cycle, verbose)
        data_to_return += [cycle_data]  #added 16J25

        if ax is not None:
            V = cycle_data[V_str]
            J = cycle_data[J_str]
    
            if colors:
                color = colors[n]
            else: color = 'k'
            
            ax.plot(V, J, color) 
            
    if ax is not None:
        ax.set_xlabel(V_str)
        ax.set_ylabel(J_str)      
    
        if saveit:
            if title == 'default':
                title == CV_data['title'] + '.png'
            fig1.savefig(title)
            
    if verbose:
        print('\nfunction \'plot_CV_cycles\' finished!\n\n')  
    
    if ax is not None:
        return data_to_return, ax
    return data_to_return


if __name__ == '__main__':
    from Data_Importing import import_data

    plt.close('all')    
    
    import_raw_data = 0
    if import_raw_data:
        EC_directory = ('/home/soren/Dropbox (Soren Scott Inc)' +
            '/Sniffer_Experiments/03_Pt_Sputtered/Data/16I28_for_Hawaii/')
        EC_file =  '03_HER_OER_C01.mpt'
        
        EC_data_0 = import_data(EC_directory + EC_file, data_type='EC')
    
    sync_metadata(EC_data_0, RE_vs_RHE=0.535, A_el=0.196)
    
    EC_data_1 = select_cycles(EC_data_0, [1,2], tzero='start')
    [ax1, ax2] = plot_vs_time(EC_data_1, colors=['r','k'])
    
    EC_data_2 = select_cycles(EC_data_0, [3,4], tzero='start')
    plot_vs_time(EC_data_2, colors=['r--','k--'], axes=[ax1, ax2])
    
    
    