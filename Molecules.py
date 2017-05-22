# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 22:36:15 2016

@author: scott

This module will define the Molecule class used to organize, access, and 
manipulate information about molecules. Objects of this class will generate and 
utilize text files in EC_MS/data/

"""

from __future__ import print_function, division
import os
import re

if os.path.split(os.getcwd())[1] == 'EC_MS':   
                                #then we're running from inside the package
    data_directory = os.getcwd() + os.sep + 'data' + os.sep
else:
    data_directory = os.getcwd() + os.sep + 'EC_MS' + os.sep + 'data' + os.sep

#for python2:
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

class Molecule:
    '''
    This class will store physical and thermodynamic data about molecules that
    we are interested in for EC_MS, as well as mass spectra and calibration 
    results for quantification. Objects of this class will link to data stored 
    in a text file in EC_MS/data/
    '''
    def __init__(self, name):
        self.attr_status = {'D':0,'kH':0,'n_el':0}
        # 0 for undefined, 1 for loaded from file, 2 for set by function
        self.name = name
        #try:
        file_name = name + '.txt'
        cwd = os.getcwd()
        os.chdir(data_directory)
        try: 
            with open(file_name,'r') as f:
                self.file_lines = f.readlines()
                if len(self.file_lines) == 0:
                    print('The file for ' + name + ' is empty! Writing new.')
                    raise FileNotFoundError
            self.reset()
        except FileNotFoundError: # I don't know the name of the error, so I'll have to actually try it first.
            print('no file found for ' + name + '. Writing new.')            
            with open(file_name,'w') as f:
                self.write_new(f)
        os.chdir(cwd)
    
    def reset(self, verbose=1):
        print('loading attributes for this ' + self.name + ' molecule fresh from file.')

        for (i,line) in enumerate(self.file_lines):
            try:
                attr = re.search(r'^\w+\s', line).group()[:-1]
                value = eval(re.search(r'\s[-]?\d+[\.]?\d*(e[-]?\d+)?\s', line).group()[1:-1])
                # I don't understand why the () in the re above can't be []
                print('got ' + attr)
                self.attr_status[attr] = 1
                setattr(self, attr, value)
            except AttributeError:    #Again, I'll have to actually try this.
                print('No data found on line ' + str(i))
                
    def write_new(self, f):
        for attr in self.attr_status:
            string = input('Enter ' + attr + ' for ' + self.name + ' as a float or anything else to skip.\n')
            try: 
                value = float(string)
                self.attr_status = 2
                setattr(self, attr, value) 
                f.write(attr + '\t=\t' + str(value) + '\n')
            except ValueError:
                print('skipped that for now.')
    
    def set_temperature(self, T):
        print('The set_temperature function is not implemented yet.')



if __name__ == '__main__':
    os.chdir('/home/soren/Dropbox (Soren Scott Inc)/Sniffer_Experiments/03_Pt_Sputtered/Analysis/16J14_time_response')
     
    if os.path.split(os.getcwd())[1] == 'EC_MS':   
                                #then we're running from inside the package
        data_directory = os.getcwd() + os.sep + 'data' + os.sep
    else:
        data_directory = os.getcwd() + os.sep + 'EC_MS' + os.sep + 'data' + os.sep    
    
    mol = Molecule('H2')    
    