

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



def select_cycles(EC_Data_0, cycles=1, verbose=1):
    ''' 
    This function selects one or more cycles from EC_Data.
    Use this before synchronizing!
    Works for both CA and CV
    '''
    if verbose:
        print('\nSelecting cycles ' + str(cycles) + ' from \'' + EC_Data_0['title'] + '\'\n')
    
    EC_Data = EC_Data_0.copy()
    
    #it looks like I actually want Ns for CA's
    if 'Ns' in EC_Data['data_cols']:
        cycle_numbers = EC_Data['Ns']
    else:
        cycle_numbers = EC_Data['cycle number']
    
    N = len(cycle_numbers)
    if type(cycles)==int:
        cycles = [cycles]
    I_keep = np.array([I for I in range(N) if cycle_numbers[I] in cycles])
    #list comprehension is awesome.

    for col in EC_Data['data_cols']:
        EC_Data[col] = EC_Data[col][I_keep]
        
    t0 = timestamp_to_seconds(EC_Data['timestamp'])
    EC_Data['tspan'] = [min(EC_Data['time/s']) + t0, max(EC_Data['time/s']) + t0]
    EC_Data['data_type'] += ' selected'   
    
    return EC_Data


def plot_CV_cycles(CV_Data, cycles=[0], Ref_vs_RHE=0.928, ax1='new',
                   saveit=1, title='default', leg=0, verbose=1, A_el=0, Colors=0):
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
        
    if Ref_vs_RHE==0 and 'Ref_vs_RHE' in CV_Data:
        Ref_vs_RHE = CV_Data['Ref_vs_RHE']
    if A_el==0 and 'A_el' in CV_Data:
        A_el = CV_Data['A_el']
    
    for n,cycle in enumerate(cycles):
        Cycle_Data = select_cycles(CV_Data, cycle, verbose)
    
        V = Cycle_Data['Ewe/V']    
        if Colors:
            color = Colors[n]
        else: color = 'k'       
           
        if Ref_vs_RHE:
            V_string = 'E vs RHE / V' 
            V = V + Ref_vs_RHE
        else:
            V_string = 'E vs ref / V'
    
        J = Cycle_Data['<I>/mA']
        if A_el:
            J = J/A_el
    
        ax1.plot(V,J,color)    
        
    if A_el:
        J_str = 'I /[mA/cm^2]'  
    else:
        J_str = 'I /[mA]'
    ax1.set_ylabel(J_str)
    
    if Ref_vs_RHE:
        V_str = 'E vs RHE / V' 
    else:
        V_str = 'E vs ref / V'        
    ax1.set_xlabel(V_str)
        
    if saveit:
        if title == 'default':
            title == EC_and_MS['title'] + '.png'
        fig1.savefig(title)
    if verbose:
        print('\nfunction \'plot_CV_cycles\' finished!\n\n')  
         

if __name__ == '__main__':
    from Data_Importing import import_data
    
    
    
    
    