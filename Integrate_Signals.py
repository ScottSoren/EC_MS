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
    from Plotting import plot_masses
else:                           #then we use relative import
    from .Plotting import plot_masses


def integrate_transient(x, y, tspan=None, t_transient=None, t_steady='half',
                        ax=None, title=None, colors=['r', 'b', 'g'], 
                        verbose=True):
    '''
    This will return seperate values for the transients and steady-states of a
    a certain compound, based on extrapolating the average signal after t_transient
    to the interval before t_transient and subtracting.
    '''
    if ax == 'new':
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    if tspan is None:
        tspan = [x[0], x[-1]]
    '''
    if t_transient is None:
        t_transient = tspan
    elif t_transient == 'half':
        t_transient = [tspan[0], (tspan[0] + tspan[-1])/2]
    '''
    if t_steady == 'half':
        t_steady = [(tspan[0] + tspan[-1])/2, tspan[1]]
        
    I_int = [I for (I, x_I) in enumerate(x) if tspan[0]<=x_I<=tspan[-1]]
#    I_transient = [I for (I, x_I) in enumerate(x) if t_transient[0]<=x_I<=t_transient[-1]]
    I_steady = [I for (I, x_I) in enumerate(x) if t_steady[0]<x_I<=t_steady[-1]]

    x_int = x[I_int]
    y_int = y[I_int]
#    x_transient = x[I_transient]
#    y_transient = y[I_steady]
    x_steady = x[I_steady]
    y_steady = y[I_steady]
    
    base = np.trapz(y_steady, x_steady) / (x_steady[-1] - x_steady[0])
 
    y_zero = np.zeros(np.shape(x_int))   
    y_base = y_zero + base    
    y_s = np.minimum(y_int, base)
    y_t = np.maximum(y_int-base, 0)
    
    steady = np.trapz(y_s, x_int)
    transient = np.trapz(y_t, x_int)
    
    if ax is not None:
        ax.fill_between(x_int, y_int, y_zero, where=y_int>y_zero,
                facecolor=colors[1], interpolate=True)
        ax.fill_between(x_int, y_int, y_base, where=y_int>y_base,
                facecolor=colors[2], interpolate=True)
        ax.plot(x, y, color=colors[0])
        if title is not None:
            ax.set_title(title)
            
    if verbose:
        if title is not None:
            print(title)
        print('\tsteady = ' + str(steady) + '\n\ttransient = ' + str(transient))
    
    return steady, transient

if __name__ == '__main__':
    pass
    