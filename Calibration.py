# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:50:10 2016
Most recently edited: 16J27

@author: scott




"""


from __future__ import division, print_function

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize, curve_fit
from scipy.integrate import odeint

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    import Chem
    from Molecules import Molecule
else:                           #then we use relative import
    from . import Chem
    from .Molecules import Molecule


   
   
if __name__ == '__main__':

    
    from Data_Importing import import_data
    from Combining import synchronize
    from EC import select_cycles, plot_CV_cycles, plot_vs_time
    from Plotting import plot_masses_and_I
    from Time_Response import stagnant_pulse, fit_step
  


        ##
    
    
    
    



