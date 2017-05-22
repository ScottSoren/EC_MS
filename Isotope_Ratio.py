# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 22:50:34 2016
last edited: 16I23

@author: Scott

This module has functions designed to analyze OER in EC-MS with isotopic labeling
"""

# make python2-compatible:
from __future__ import print_function
from __future__ import division

import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy


def predict_M34(MS_Data, subtract_background = 1):
    x_34 = 1/2 * (MS_Data['M32-x'] + MS_Data['M36-x'])
    
    y_32 = MS_Data['M32-y']
    y_36 = MS_Data['M36-y']
    
    if subtract_background:
        b_32 = min(y_32)     
        y_32 = y_32 - b_32        
        b_36 = min(y_36)
        y_36 = y_36 - b_36
    
    r = y_36/y_32
    y_34 = y_32*2*np.sqrt(r)
    if subtract_background:
        y_34 = y_34 + min(b_32,b_36)
    
    MS_Data['predicted M34-x'] = x_34
    MS_Data['predicted M34-y'] = y_34
    MS_Data['data_cols'].append('predicted M34-x') 
    MS_Data['data_cols'].append('predicted M34-y') 

    return MS_Data



if __name__ == '__main__':
    
    from Data_Importing import import_data
    from Combining import numerize, synchronize, time_cut, plot_masses, plot_masses_and_I    

    plt.close()    
    
    default_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))  

    #MS_Data = numerize(import_data(MS_file,data_type = 'MS'))
    import_raw_data = 1
    if import_raw_data:
        MS_Data_0 = import_data(default_directory, data_type='MS')
        CA_Data_0 = import_data(default_directory, data_type='EC')
        
    CV_Data = numerize(MS_Data_0)
    CA_Data = numerize(CA_Data_0)    
    
    
    CA_and_MS = synchronize([MS_Data, CA_Data], cutit = 0)    
    tspan = CA_and_MS['tspan_2']   
    tspan[1] = tspan[1] + 100 + 0*(tspan[1]-tspan[0])
    tspan[0] = tspan[0] - 20 - 0*(tspan[1]-tspan[0])  

    MS_Data_1 = time_cut(deepcopy(CA_and_MS),tspan)

    MS_Data_1 = predict_M34(MS_Data_1)    
    
    plot_masses_and_I(MS_Data_1, logplot = [1,1], overlay = 0,
                Colors = {'M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m', 'predicted M34':'--k'})
     #            Colors = {'M2':'k','M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m'})

    
    
    
    
    