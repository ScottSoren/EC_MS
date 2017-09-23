#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 00:21:47 2017

@author: scott


17G30:
There may be something time-saving here:
/home/scott/Dropbox/other_DTU/MSc/Scripts/PYTHON3/Masters Project/15L29 Report Intro
... but I can probably do better from scratch

or, better, build off someone elses work. 
how about this?
https://pypi.python.org/pypi/CoolProp/2.2.3
#or a function that looks up from hbcp or something...
"""

import numpy as np
R = 8.3144598            #gas constant / (J/(mol*K))  #I'd like to import from PhysCon, but relative import is a pain.


S0 = { # standard entropy / [J/(mol*K)]
      'H2O(g)':188.72, 'H2O(l)':69.940, 
      'ethanol(g)':282, 'ethanol(l)':161, 'ethanol(aq)':None
      }
dfH0 = { # standard enthalpies of formation / [kJ/mol]
        'H2O(g)': -241.82 ,'H2O(l)': -285.8, 
        'ethanol(g)':-235.3, 'ethanol(l)':-277
        }
dfG0 = { # standard enthalpies of formation / [kJ/mol] 
         # to be populated later
        } 

dsH0 = {'ethanol':-19.5,  #temperarily taken value for CO2 below to check that the function runs properly! Couldn't find it for ethanol
        'CO2':-19.5,  #Carroll1991
        } #solvation enthalpy at 25C

kH = {'ethanol':None}

def p_vap(mol='H2O', T=298.15, unit='Pa'):
    dH = (dfH0[mol+'(g)'] - dfH0[mol+'(l)']) * 1e3 
    dS = S0[mol+'(g)'] - S0[mol+'(l)']
    
    if unit == 'Pa':
        p0 = 1e5
    elif unit == 'bar':
        p0 = 1
    elif unit == 'mbar':
        p0 = 1000
    
    p = p0 * np.exp( -dH/(R*T) + dS/R )
    
    return p



def kH_of_T(T, kH_0=None, dsH=None, mol=None, T0=298.15):
    
    if kH_0 is None:
        kH_0 = kH[mol]
    
    if dsH is None:
        dsH = dsH0[mol]
    
    kH_T = kH_0 * np.exp( - dsH/R * (1/T0 - 1/T))
    
    return kH_T



if __name__ == '__main__':
    print(p_vap('H2O', T=298.15, unit='mbar'))
    
    
    
    