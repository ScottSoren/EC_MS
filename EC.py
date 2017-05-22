

# -*- coding: utf-8 -*-
"""
created 16I15
EC.py
last edited 16I23

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
    from Combining import timestamp_to_seconds
else:                           #then we use relative import
    from .Combining import timestamp_to_seconds



def select_cycles(EC_data_0, cycles=1, verbose=1):
    ''' 
    This function selects one or more cycles from EC_data_0.
    Use this before synchronizing!
    Works for both CA and CV
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
        EC_data[col] = EC_data[col][I_keep]
        
    t0 = timestamp_to_seconds(EC_data['timestamp'])
    EC_data['tspan'] = [min(EC_data['time/s']) + t0, max(EC_data['time/s']) + t0]
    EC_data['data_type'] += ' selected'   
    
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

def plot_CV_cycles(CV_data, cycles=[0], RE_vs_RHE=0, ax1='new',
                   saveit=0, title='default', leg=0, verbose=1, A_el=0, colors=0):
    '''
    plots a subset of cycles in a CV
    '''
    if verbose:
        print('\n\nfunction \'plot_CV_cycles\' at your service!')
    
    if ax1 == 'new':
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111) 
    
    if type(cycles)==int:
        cycles = [cycles]  
        
    if RE_vs_RHE==0 and 'RE_vs_RHE' in CV_data:
        RE_vs_RHE = CV_data['RE_vs_RHE']
    if A_el==0 and 'A_el' in CV_data:
        A_el = CV_data['A_el']
    
    for n,cycle in enumerate(cycles):
        cycle_data = select_cycles(CV_data, cycle, verbose)
    
        V = cycle_data['Ewe/V']    
        if colors:
            color = colors[n]
        else: color = 'k'
        
        J = cycle_data['<I>/mA']
        if A_el:
            J = J/A_el
    
        ax1.plot(V,J,color)    
        
    if A_el:
        J_str = 'I /[mA/cm^2]'  
    else:
        J_str = 'I /[mA]'
    ax1.set_ylabel(J_str)
    
    if RE_vs_RHE:
        V_str = 'E vs RHE / V' 
    else:
        V_str = 'E vs ref / V'        
    ax1.set_xlabel(V_str)
        
    if saveit:
        if title == 'default':
            title == CV_data['title'] + '.png'
        fig1.savefig(title)
    if verbose:
        print('\nfunction \'plot_CV_cycles\' finished!\n\n')  
         

if __name__ == '__main__':
    from Data_Importing import import_data
    
    
    
    
    