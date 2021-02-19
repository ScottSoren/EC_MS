# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 11:39:18 2021

@author: scott
"""
import json, re
from .Molecules import Molecule

def mdict_from_SI2020_calibration_file(calibration_file):
    """Return a molecule dict from an old Spectro Inlets calibration file"""
    mdict = {}
    
    with open(calibration_file, "r") as f:
        calibration_dict = json.load(f)
    
    for name in calibration_dict["mol_list"]:
        real_name = calibration_dict["real_names"].get(name, name)
        m = Molecule(real_name)
        m.F_mat = {}
        
        primary_match = re.search(r"M[0-9]+", name)
        if primary_match:
            m.primary = primary_match.group(0)
        
        if name in calibration_dict["F"]:
            for mass, F_M in calibration_dict["F"][name].items():
                if F_M == 0:
                    continue  # EC_MS can't handle sensitivity=0
                m.F_mat[mass] = F_M
                if mass == m.primary:
                    m.F = F_M
                    m.F_cal = F_M
                
        mdict[name] = m
    
    return mdict
    