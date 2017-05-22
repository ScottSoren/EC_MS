# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 22:50:34 2016
Most recently edited: 16I23

@author: Scott

This module has functions for integrating and averaging signals over specified
time frames or cycles, mostly for pulsed electrolysis experiments.
"""
# make python2-compatible:
from __future__ import print_function
from __future__ import division

import numpy as np
from matplotlib import pyplot as plt
import os

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from Combining import plot_masses
else:                           #then we use relative import
    from .Combining import plot_masses



def get_cumulative(Data, data_cols,  tspan = 0, verbose = 0, value_type = 'integral'):
    '''gets an value integrated or averaged over tspan for a variable in a dataset'''
    if type(data_cols) == str:
        data_cols = [data_cols]
    values = []
    for data_col in data_cols:
        if verbose:
            print('\nGetting ' + value_type + ' for ' + data_col + ' in ' + Data['title']) 
        if data_col[0] == 'M':         #dealing with QMS data          
            xcol = data_col + '-x'
            ycol = data_col + '-y'
        else:                       #dealing with EC data
            xcol = 'time/s'
            ycol = data_col
        x0 = Data[xcol]
        y0 = Data[ycol]    
        if tspan:
            if len(tspan)==1: #then it is the length of time to integrate from start
                tspan = [min(x0), min(x0) + tspan]
            I_keep = np.array([I for I in range(len(x0)) if x0[I]>=tspan[0] if x0[I]<=tspan[1]])
            x = x0[I_keep]      #x.copy() here was much slower.
            y = y0[I_keep]  
        if value_type == 'integral':    
            value = np.trapz(y,x)
        elif value_type == 'average':
            value = np.trapz(y,x)/(x[-1] - x[0])   
        if verbose:
            print(value_type + ' = ' + str(value) + '\n')
        values += [value]
    if len(values) == 1:
        return values[0]
    return values


def step_vs_potential(CA_and_MS, data_cols, verbose=1, value_type_1='integral',
                       cycles='all', excluded_cycles='none', t_int=-1, extra=0, plotit=0):
    '''
    returns the cumulative values for a given set of variables over a set of CA
    steps
    '''
    
    title = CA_and_MS['title']
    if verbose:
        print('\nGetting integrated ' + str(data_cols) + ' vs potential for cycles in ' + title)
    
    if cycles == 'all':
        N_cycles = CA_and_MS['Ns'][-1] + 1 #since CA's start at cycle zero
        cycles = range(int(N_cycles))
    if type(excluded_cycles) == list:
        for excluded_cycle in excluded_cycles:
            cycles.remove(excluded_cycle)
            
    N_cycles = len(cycles)
    if type(data_cols) is not list:
        data_cols = [data_cols]
    N_cols = len(data_cols)
    Potentials = np.zeros(N_cycles)
    Values = np.zeros([N_cycles,N_cols])
    tspans = np.zeros(N_cycles)
    
   
    if type(t_int) == float or type(t_int) == int:    #t_int is integration time, for specified for each step or once for all.
        if t_int == -1:     #in this case all steps will be integrated fully
            t_int = -1*np.ones([N_cycles])        
        else: #in this case all steps will be integrated for the same length of time
            t_int = t_int*np.ones([N_cycles])
    
    cycle_numbers = CA_and_MS['Ns']
    t = CA_and_MS['time/s']
    for n_cycle, c in enumerate(cycles):
        t_cycle = np.array([t[I] for I in range(len(cycle_numbers)) if cycle_numbers[I]==c])
        #I love list comprehension
        tspan0 = [t_cycle[0], t_cycle[-1]]
        
        tspan = tspan0.copy()
        if t_int[n_cycle]>=0:
            if extra:
                tspan[1] = tspan[1]+t_int[n_cycle]   
            else:
                tspan[1] = tspan[0]+t_int[n_cycle] 
        tspans[n_cycle] = t[np.where(t<tspan[1])[0][-1]] - t[np.where(t>tspan[0])[0][0]] #gives the actual timespan,taking into account the QMS's poor resolution
                
            #so that I can specify to only integrate over the first t_int seconds of the rest period
        potential = get_cumulative(CA_and_MS, 'Ewe/V', tspan0, value_type = 'average') #average potential only over the pulse (tspan0)
        value = get_cumulative(CA_and_MS, data_cols, tspan, value_type = value_type_1) #integrate signal over pulse plus tail (tspan1)
        #value of each integral at that potential
        if verbose:
            print('cycle ' + str(c) + ' potential ' + str(potential) + ' value ' + str(value))
            if plotit:
                plot_masses(CA_and_MS, tspan=tspan, logplot=1, verbose=1,
                            Colors = {'M2':'b','M15':'r','M26':'g','M28':'0.5','M32':'k'}, 
                            ax1='new', saveit=0, leg=0)
        Potentials[n_cycle] = potential
        Values[n_cycle,:] = value
        
    if verbose:
        print('\tgot ' + str(len(Potentials)) + ' data points!\n' )
    

        
    return (Potentials, Values, tspans)


def integrate_transient(CA_and_MS, mass='M15', t_int=15, t_transient=20, cycles='all', verbose=1, plotit=0):
    '''
    This will return seperate values for the transients and steady-states of a
    a certain compound, based on extrapolating the average signal after t_transient
    to the interval before t_transient and subtracting.
    '''
    (V, s_total, tspan) = step_vs_potential(CA_and_MS, mass, cycles=cycles, t_int=0, extra=1, plotit=plotit) 
    s_total = s_total[:,0]      #I still really don't like numpy arrays. Makes me miss Matlab.
    (_, s_early, tspan_early) = step_vs_potential(CA_and_MS, mass, cycles=cycles, t_int=t_transient, extra=0, plotit=plotit) 
    s_early = s_early[:,0]
    s_late = s_total - s_early 
    s_steady = s_late * tspan / (tspan - tspan_early) 
    s_transient = s_total - s_steady 

    
    return(s_steady, s_transient)
    


def integrate_pulses(CA_and_MS, masses, verbose=1,
                     cycles='all', excluded_cycles='none', t_tail=15):
    '''
    Returns a set of integrated QMS signals vs potential for Scott and Daniel's
    pulsed CORR experiments on the sniffer... a bit ugly.
    '''
    title = CA_and_MS['title']
    if verbose:
        print('\n\nfunction \'integrate pulses\' at your service!\n' +
        'Integrating pulses on \'' + title + '\'')
    data_cols = ['I/mA'] + masses

    N_cycles = CA_and_MS['Ns'][-1] + 1 #since cycle_number starts at 0
    t_int = -1* np.ones([N_cycles]) #pulses t
    t_int[::2] = t_tail #this ensures that the rest periods are only integrated for the first few seconds
    
    (E, A, tspans) = step_vs_potential(CA_and_MS, data_cols, cycles=cycles, excluded_cycles=excluded_cycles, t_int=t_int)
   
   #first cycle is rest, pulse is at second cycle (index 0). Total number of pulses is odd
    dQ = A[1::2,0]                #only keep charge from pulses
    E_pulse = E[1::2]                #only keep potential from pulses
    t_pulses = tspans[1::2]
    A_M = A[1::2,1:] + A[2::2,1:] #add tail to pulse signal for products
    results_title = title + '_results'
    Pulse_Results = {'title':results_title, 'masses': masses,
                    'E vs ref':E_pulse, 'dQ': dQ, 't_pulses': t_pulses, 'tspans': tspans}
                    
    if 'Ref_vs_RHE' in CA_and_MS:
        Pulse_Results['E vs RHE'] = Pulse_Results['E vs ref'] + CA_and_MS['Ref_vs_RHE']
    for i, mass in enumerate(masses):
        Pulse_Results[mass] = A_M[:,i]
    
    if verbose:
        print('function \'integrate pulses\' finished!\n\n')
    return Pulse_Results


def plot_pulse_integrals(Pulse_Results, Colors):    

    x = Pulse_Results['E vs RHE']    
    
    fig1 = plt.figure()
    
    ax1 = fig1.add_subplot(211)
    ax1.plot(x,-1*Pulse_Results['dQ'],'k.',label = 'charge')
    ax1.set_yscale('log')
    ax1.set_xlabel('E vs RHE / V')
    ax1.set_ylabel('integrated signal')
    
    ax2 = fig1.add_subplot(212)
    for mass in Pulse_Results['masses']:    
        spec = Colors[mass] + '.'
        y = Pulse_Results[mass]
        ax2.plot(x, y, spec, label=mass)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.set_ylabel('dQ / C')    
    
    return (ax1, ax2)


if __name__ == '__main__':
    import os
    from Data_Importing import import_data
    from Combining import numerize, synchronize, plot_masses
    #plt.close()    
    
    default_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))  

    #MS_Data = numerize(import_data(MS_file,data_type = 'MS'))
    import_raw_data = 0
    if import_raw_data:
        MS_Data_0 = import_data(default_directory, data_type='MS')
        CA_Data_0 = import_data(default_directory, data_type='EC')
    
    MS_Data = numerize(MS_Data_0)
    CA_Data = numerize(CA_Data_0)  
    
    CA_and_MS = synchronize([MS_Data, CA_Data], cutit = 1)    
    tspan = CA_and_MS['tspan_2']   
    #tspan[1] = tspan[1] + 100 + 0*(tspan[1]-tspan[0])
    #tspan[0] = tspan[0] - 20 - 0*(tspan[1]-tspan[0])     
    
    Colors = {'M2':'b', 'M15':'r', 'M26':'g'}
    masses = ['M2', 'M15', 'M26']

    Pulse_Results = integrate_pulses(CA_and_MS, masses, verbose = 1)   

    ref_vs_RHE = 0.918
    Pulse_Results['E vs RHE'] = Pulse_Results['E vs ref'] + ref_vs_RHE    
    
    (ax1,ax2) = plot_pulse_integrals(Pulse_Results, Colors)
    