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

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    import Chem
    from EC import plot_CV_cycles, CV_difference, sync_metadata
    from Molecules import Molecule
else:                           #then we use relative import
    from . import Chem
    from .EC import plot_CV_cycles, CV_difference, sync_metadata
    from .Molecules import Molecule



def ML_strip_cal(CV_and_MS, cycles=[1,2], t_int=200,
             mol='CO2', mass='primary', n_el=None, 
             Vspan=[0.5, 1.0], redox=1, 
             ax='two', title='default', verbose=True):
    '''
    Determines F_cal = Q_QMS / n_electrode by integrating a QMS signal over
    tspan, assuming the starting value is background; and
    integrating over vspan the difference between two CV cycles and converting 
    that to a molar amount. 
    Returns a partially populated calibration dictionary. 
    The calibration factor is calibration['F_cal']
    '''
    if verbose:
        print('\n\ncalibration function \'ML_strip_cal\' at your service!\n')
        
    if ax == 'two':
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(211)
        ax2 = fig1.add_subplot(212)
    elif ax is list:
        ax1 = ax[0]
        ax2 = ax[1]
    elif ax is None:
        ax1 = None
        ax2 = None
    else:
        ax1 = ax
        ax2 = None
        
    if type(mol) is str:
        mol = Molecule(mol, writenew=False)
    name = mol.name
    if n_el is None:
        n_el = mol.n_el
    if mass == 'primary':
        mass = mol.primary
    if np.size(t_int) == 1:
        t_int = CV_and_MS['tspan_2'][0] + np.array([0, t_int])
    if title == 'default':
        title = name + '_' + mass
    
    cycles_data, ax1 = plot_CV_cycles(CV_and_MS, cycles, ax=ax1, title=title)
    Q_diff, diff = CV_difference(cycles_data, Vspan=Vspan, redox=redox, ax=ax1)
    ax1.set_title(title)    
    
    n_mol = Q_diff / (Chem.Far * n_el)   
    t = diff[0][0]
    J_diff = diff[2]
    A_el = CV_and_MS['A_el']    
    
    #Q_diff seemed to be having a problem, but turned out to just be because 
     #I forgot to remove_delay(). 
    # Now everything checks out, as can be seen here:
    '''
    Q_diff1 = A_el * np.trapz(J_diff, t) * 1e-3   #factor 1e-3 converts mA to A
    print('Q_diff = ' + str(Q_diff) +'\nQ_diff1 = ' + str(Q_diff1) + 
          '\nratio = ' + str(Q_diff1/Q_diff))
    '''
    x = CV_and_MS[mass + '-x']
    y = CV_and_MS[mass + '-y']
    I_keep = [I for (I, x_I) in enumerate(x) if t_int[0] < x_I < t_int[1]]
    x = x[I_keep]
    y = y[I_keep]
    
    
    background = min(y)  #very simple background subtraction
    Q_QMS = np.trapz(y - background, x)
    F_cal = Q_QMS / n_mol
    
    y_el = J_diff * A_el/(Chem.Far * n_el) * F_cal * 1e-3 
        #factor 1e-3 converts mA to A
    
    if ax2 is not None:
        ax2.plot(x,y, 'k-')
        ax2.fill_between(x, background*1e9, y*1e9, where=y>background, 
                         facecolor='g', interpolate=True)
        ax2.plot(t, y_el*1e9, 'r--')    #Without mass transport     
        ax2.set_xlabel('time / s')
        ax2.set_ylabel('signal / nA')
#        ax2.set_yscale('log')
        
    
    print(('QMS measured {0:5.2e} C of charge at M44 for {1:5.2e} mol ' + name + '.\n' + 
            'Calibration factor for CO2 at M44 is {2:5.2e} C / mol.'
            ).format(Q_QMS, n_mol, F_cal))    
    
    calibration = {'type': 'ML_strip'}
    calibration['raw_data'] = CV_and_MS['title']    
    calibration['mass'] = mass
    calibration['n_mol'] = n_mol
    calibration['Q_el'] = Q_diff
    calibration['Q_QMS'] = Q_QMS    
    calibration['F_cal'] = F_cal
    
    if verbose:
        print('\ncalibration function \'ML_strip_cal\' finished!\n\n')    

    if ax is None:
        return calibration
    if ax2 is None:
        ax = ax1
    else:
        ax = [ax1, ax2]
    return calibration, ax



def steady_state_cal(CA_and_MS, t_int='half',
                     mol='CO2', mass='primary', n_el=None, 
                     ax='new', title='default', verbose=True,
                     background='min'):
    if verbose:
        print('\n\ncalibration function \'steady_state_cal\' at your service!\n')
    if type(mol) is str:
        mol = Molecule(mol, writenew=False)
    name = mol.name
    if n_el is None:
        n_el = mol.n_el
    if mass == 'primary':
        mass = mol.primary
    if t_int == 'half':
        t_int = (CA_and_MS['tspan_2'][1] + np.array(CA_and_MS['tspan_2'])) / 2
    elif t_int == 'all':
        t_int =  np.array(CA_and_MS['tspan_2'])
    elif np.size(t_int) == 1:
        t_int = CA_and_MS['tspan_2'][1] + np.array([-t_int, 0])    
        #by default integrate for time t_int up to end of interval         
    if title == 'default':
        title = name + '_' + mass
    x = CA_and_MS[mass + '-x']
    y = CA_and_MS[mass + '-y']
    if background == 'min':
        background = min(y)
    elif background is None:
        background = 0
    I_keep = [I for (I, x_I) in enumerate(x) if t_int[0] < x_I < t_int[1]]
    x_r = x[I_keep]
    y_r = y[I_keep]
    Q_QMS = np.trapz(y_r - background, x_r)          #integrated signal in C

    V_str, J_str = sync_metadata(CA_and_MS)
    t = CA_and_MS['time/s']
    J = CA_and_MS[J_str]
    A_el = CA_and_MS['A_el']    
    
    I_keep = [I for (I, t_I) in enumerate(t) if t_int[0] < t_I < t_int[1]]
    t_r = t[I_keep]
    J_r = J[I_keep]
    Q_el = A_el * np.trapz(J_r, t_r) * 1e-3 # total electrode charge passed in C
    n_mol = Q_el / (Chem.Far * n_el)
    
    F_cal = Q_QMS / n_mol
    
    y_el = J * A_el/(Chem.Far * n_el) * F_cal * 1e-3 
    # expected QMS signal without mass transport etc
    
    print(('QMS measured {0:5.2e} C of charge at M44 for {1:5.2e} mol ' + name + '.\n' + 
            'Calibration factor for CO2 at M44 is {2:5.2e} C / mol.'
            ).format(Q_QMS, n_mol, F_cal))     

    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    if ax is not None:
        ax.plot(x, y, 'k-')
        ax.plot(t, y_el + background, 'r--')
        ax.set_title(title)
        
    calibration = {'type': 'steady_state'}
    calibration['raw_data'] = CA_and_MS['title']    
    calibration['mass'] = mass
    calibration['n_mol'] = n_mol
    calibration['Q_el'] = Q_el
    calibration['Q_QMS'] = Q_QMS    
    calibration['F_cal'] = F_cal         

    if verbose:
        print('\ncalibration function \'steady_state_cal\' finished!\n\n')  
    
    return calibration



def LCA_cal():
    calibration = {'type': 'LCA'}
    return calibration
   
   
if __name__ == '__main__':
    '''
    
    from Data_Importing import import_data
    from Combining import synchronize
    from EC import select_cycles, plot_CV_cycles, plot_vs_time
    from Plotting import plot_masses_and_I
    from Time_Response import stagnant_pulse, fit_step
  


        ##
    
    
    
    '''   



