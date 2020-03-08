# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:28:25 2020

@author: scott
"""

import os, re

data_directory = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data"


def get_NIST_spectrum(mol):
    """
    a parser for NIST-exported .jdx files
    """
    data_folder = data_directory + os.sep + "NIST_spectra_data"
    if type(mol) is not str:
        try:
            mol = mol.real_name
        except AttributeError:
            mol = mol.name

    file_list = os.listdir(data_folder)

    try:
        file = next(f for f in file_list if re.search("^" + mol, f))
    except StopIteration:
        print("WARNING!!! No Spectrum available for " + mol)
        raise FileNotFoundError
    # ^ file-extension-ambiguous because I might forget and save them as .txt at some point

    with open(data_folder + os.sep + file) as f:
        lines = f.readlines()

    in_data = False
    spectrum = {}
    for line in lines:
        if "END" in line:
            break
        if in_data:
            # print(line) # debugging
            mass_values = line.strip().split(" ")
            for mass_value in mass_values:
                mass, value = mass_value.split(",")
                mass = "M" + mass
                value = eval(value)
                spectrum[mass] = value
        elif "PEAK TABLE" in line:
            in_data = True

    return spectrum
