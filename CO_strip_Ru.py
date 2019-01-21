#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 18:41:11 2018

@author: scott
"""
import os
import pickle
from matplotlib import pyplot as plt

from EC_MS import sync_metadata, make_selector, select_cycles, plot_experiment
from EC_MS import get_shunt_current_line, correct_shunt, correct_ohmic_drop
from EC_MS import plot_vs_potential

plt.close('all')

data_directory = os.path.expanduser('~/Dropbox/DTU-MIT RuO2/Data/ECMS/')
A_el = 0.196

measurements = [
        ('Melih1', '18I13_Reshma1D_18O', '01', '18O', 0.720, 500),
        ()
                ]




for name, folder, tag, isotope, RE_vs_RHE, R in measurements:
    
    F = data_directory + folder
    pkl_file_name = [f for f in  os.listdir(F) if f[0:2]==tag and f[-4:]=='.pkl'][0]
    with open(F + os.sep + pkl_file_name, 'rb') as pkl:
        data = pickle.load(pkl)
        last_pkl_file_name = pkl_file_name
    
    # calibrate data
    #t_str = trigger_cal(data) # ligns up EC and MS time axes
    t_str = 'time/s'
    V_str, J_str = sync_metadata(data, RE_vs_RHE=RE_vs_RHE, A_el=A_el) 
        # ^calibrates EC data
    sel_str = make_selector(data)       
    
    cycle = select_cycles(data, cycles=[0,1,2,3], cycle_str=sel_str)
    if True:
        ax0 = plot_experiment(cycle, masses=[['M36','M34','M32'],['M44','M46','M48']], logplot=False)
        plt.savefig('CO2_' + name + '_first_cycles.png')
        #ax = plot_vs_potential(cycle, masses=[['M36','M34','M32'],['M44','M46','M48']], logplot=False)
        #plt.savefig('CO2_' + name + '_first_cycles_CV.png')
    else:
        ax = [None,None]
    pfit, t_f = get_shunt_current_line(data, V_DL=[1.21, 1.31], t_i=150, ax=ax[1])
    #correct_shunt(data, pfit=pfit)
    V_str = correct_ohmic_drop(data, R_ohm=R)
    
    V_str, J_str = sync_metadata(data)
    
    cycles_data = select_cycles(data, cycles=[1,2,3])