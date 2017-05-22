#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 16:13:26 2017

To have all the tools used for the Errorbars script in versitile, 
easily accessible form.

@author: scott
"""

from matplotlib import pyplot as plt
import numpy as np




def fill_with(quantitydict, value):
    emptydict ={}
    for (key, val) in quantitydict.items():
#        print(str(key) + ' ' + str(value))
        if type(val) is dict:
            emptydict[key] = get_empty(val)
        else:
            if type(value) in [list, dict]:
                value = value.copy() # otherwise they get linked... 
            emptydict[key] = value
    return emptydict

def get_empty(quantitydict):
    return fill_with(quantitydict, value=[])
    
def add_datapoint(source, target, index=None):
#    print(str(source))
    for key, value in source.items():

        if type(value) is dict:
            add_datapoint(value, target[key], index)
        elif key in target.keys():
            
            if index is None:
                target[key] += [value]
            else:
                target[key] += [value[index]]
#                print('adding ' + str(value[index]) + ' to ' + str(key))

def values_to_stats(values):
    stats={}
    for key, value in values.items():
        #print(key)
        if type(value) is dict:
            stats[key] = values_to_stats(value)
        elif type(value) is list or type(value) is np.ndarray:
            mean = np.mean(value)
            #std = 0
            std = np.std(value)
            stats[key] = [mean, std]
    return stats
    


def get_mlu(stat): # mlu stands for for: mean, [lower, upper]
    if type(stat) is not list:
        return stat, None
    elif len(stat) == 2:
        mean = stat[0]
        upper = mean + stat[1]
        lower = mean - stat[1]
    elif len(stat) == 3:
        lower = stat[0]
        mean = stat[1]
        upper = stat[2]
    else:
        print('need stats of length 2 or 3 for errorbars')
        raise ValueError
    return mean, [lower, upper]


def plot_errorbar(xstat, ystat, ax=plt.gca(), 
                  spec='.', color='k', markersize=None):
    if markersize is None:
        if spec == '.':
            markersize = 15
        else:
            markersize = 5
    x, x_lu = get_mlu(xstat)
    y, y_lu = get_mlu(ystat)   
    ax.plot(x, y, spec, markersize=markersize, color=color)
    if x_lu is not None:
        ax.plot([x_lu[0], x_lu[1]], [y, y], '|-', color=color)
    if y_lu is not None:
        ax.plot([x, x], [y_lu[0], y_lu[1]], '_-', color=color)    

        
def plot_errorbars_y(stats, x='outer', ax='new', 
                     colors='k', Vrange=None):
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
#    print(type(stats))
    if type(stats) is not dict:
        if Vrange is None or Vrange[0] <= x <= Vrange[1]:
            plot_errorbar(x, stats, ax=ax, color=colors)
#        print('I should have just plotted something.')
        return ax
    
    if (x not in ['outer', 'inner'] and type(colors) is not dict):
        colors = fill_with(stats, colors)
    if (x not in ['outer', 'inner'] and type(Vrange) is not dict):
        colors = fill_with(stats, Vrange)

    for key, val in stats.items():
        if x=='outer':
            x_val = key
            color_val = colors
            Vrange_val = Vrange
        elif x=='inner':
            print('errorbars: x=\'inner\' not yet implemented.')
            pass
        else:
            x_val = x
            color_val = colors[key]
            Vrange_val = Vrange[key]
        plot_errorbars_y(val, x=x_val, ax=ax, colors=color_val, Vrange=Vrange_val)
        
    return ax  


