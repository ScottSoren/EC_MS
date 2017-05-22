# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 22:50:34 2016

@author: soren
"""


from Data_Importing import import_data
from EC_MS import numerize, synchronize, time_cut, plot_masses_and_I, plot_masses
import os
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy


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


def step_vs_potential(CA_and_MS, data_cols, verbose = 1, value_type_1 = 'integral',
                       cycles = 'all', excluded_cycles = 'none', t_int = -1):
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
    N_cols = len(data_cols)
    Potentials = np.zeros(N_cycles)
    Values = np.zeros([N_cycles,N_cols])
    tspans = np.zeros(N_cycles)
    
   
    if type(t_int) == float:    #t_int is integration time, for specified for each step or once for all.
        if t_int == -1:     #in this case all steps will be integrated fully
            t_int = -1*np.ones([N_cycles])        
        else: #in this case all steps will be integrated for the same length of time
            t_int = t_int*np.ones([N_cycles])
    
    cycle_numbers = CA_and_MS['Ns']
    t = CA_and_MS['time/s']
    for n_cycle, c in enumerate(cycles):
        t_cycle = np.array([t[I] for I in range(len(cycle_numbers)) if cycle_numbers[I]==c])
        #I love list comprehension
        tspan = [t_cycle[0], t_cycle[-1]]
        if t_int[n_cycle]>=0:
            tspan[1] = tspan[0]+t_int[n_cycle]      
            #so that I can specify to only integrate over the first t_int seconds of the rest period
        potential = get_cumulative(CA_and_MS, 'Ewe/V', tspan, value_type = 'average')
        value = get_cumulative(CA_and_MS, data_cols, tspan, value_type = value_type_1) 
        #value of each integral at that potential
        if verbose:
            print('cycle ' + str(c) + ' potential ' + str(potential) + ' value ' + str(value))
        Potentials[n_cycle] = potential
        Values[n_cycle,:] = value
        tspans[n_cycle] = tspan[1] - tspan[0]
    if verbose:
        print('\tgot ' + str(len(Potentials)) + ' data points!\n' )
        
    return (Potentials, Values, tspans)



def integrate_pulses(CA_and_MS, masses, verbose = 1,
                     cycles = 'all', excluded_cycles = 'none', t_tail = 15):
    '''
    Returns a set of integrated QMS signals vs potential for Scott and Daniel's
    pulsed CORR experiments on the sniffer
    '''
    title = CA_and_MS['title']
    if verbose:
        print('\n\nfunction \'integrate pulses\' at your service!\n' +
        'Integrating pulses on \'' + title + '\'')
    data_cols = ['I/mA'] + masses

    N_cycles = CA_and_MS['Ns'][-1] + 1 #since cycle_number starts at 0
    t_int = -1* np.ones([N_cycles]) #pulses t
    t_int[::2] = t_tail #this ensures that the rest periods are only integrated for the first few seconds
    
    (E, A, tspans) = step_vs_potential(CA_and_MS,data_cols, cycles=cycles, excluded_cycles=excluded_cycles, t_int=t_int)
   
   #first cycle is rest, pulse is at second cycle (index 0). Total number of pulses is odd
    dQ = A[1::2,0]                #only keep charge from pulses
    E_pulse = E[1::2]                #only keep potential from pulses
    t_pulses = tspans[1::2]
    A_M = A[1::2,1:] + A[2::2,1:] #add tail to pulse signal for products
    results_title = title + '_results'
    Pulse_Results = {'title':results_title, 'masses': masses,
                    'E vs ref':E_pulse, 'dQ': dQ, 't_pulses': t_pulses, 'tspans': tspans}
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
        ax2.plot(x, y, spec, label = mass)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.set_ylabel('dQ / C')    
    
    return (ax1, ax2)


if __name__ == '__main__':
    
    #plt.close()    
    
    directory_name = os.getcwd() + '/'    

    MS_file = directory_name + 'QMS_data_16F28'
    CA_file = directory_name + 'measurement_2_11_CA_C01.mpt'

    #MS_Data = numerize(import_data(MS_file,data_type = 'MS'))
    import_raw_data = 0
    if import_raw_data:
        MS_Data_0 = import_data(MS_file,data_type = 'MS')
        CA_Data_0 = import_data(CA_file)
    
    MS_Data = numerize(MS_Data_0)
    CA_Data = numerize(CA_Data_0)  
    
    CA_and_MS = synchronize([MS_Data, CA_Data], cutit = 1)    
    tspan = CA_and_MS['tspan_2']   
    #tspan[1] = tspan[1] + 100 + 0*(tspan[1]-tspan[0])
    #tspan[0] = tspan[0] - 20 - 0*(tspan[1]-tspan[0])     
    
    Colors = {'M2':'b', 'M15':'r', 'M27':'g', 'M31':'m'}

    #Pulse_Results = integrate_pulses(CA_and_MS, masses, verbose = 1)   

    ref_vs_RHE = 0.918
    Pulse_Results['E vs RHE'] = Pulse_Results['E vs ref'] + ref_vs_RHE    
    
    (ax1,ax2) = plot_pulse_integrals(Pulse_Results, Colors)
    