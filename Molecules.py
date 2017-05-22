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
import numpy as np
from matplotlib import pyplot as plt

if os.path.split(os.getcwd())[1] == 'EC_MS':   
                                #then we're running from inside the package
    import Chem
    from Object_Files import structure_to_lines, lines_to_dictionary, lines_to_structure, date_scott
    data_directory = os.getcwd() + os.sep + 'data' + os.sep
else:
    from . import Chem
    from .Object_Files import structure_to_lines, lines_to_dictionary, lines_to_structure, date_scott
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
    in a text file in ./data/
    '''
    def __init__(self, name, writenew=True):
        self.name = name
        self.cal = {}
        self.__str__ = '<' + name + ', instance of EC_MS class \'Molecule\'>'
        self.M = Chem.Mass(name) #molar mass from chemical formula
        self.primary = None     #the primary mass to measure at
        self.calibrations = []      #will store calibration data
        self.attr_status = {'D':0,'kH':0,'n_el':0, 'name':2}
        # 0 for undefined, 1 for loaded from file, 2 for set by function
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
            self.file_lines = ['name: ' + self.name, 'M\t=\t' + str(self.M)] + self.file_lines
        except FileNotFoundError: # I don't know the name of the error, so I'll have to actually try it first.
            print('no file found for ' + name)  
            if writenew:
                print('Writing new.\n')
                self.write(self.write_new)
        os.chdir(cwd)
    
    
    def write(self, a, *args, **kwargs):
        cwd = os.getcwd()
        os.chdir(data_directory)
        file_name = self.name + '.txt'
        with open(file_name, 'a') as f:
            if callable(a):
                a(f, *args, **kwargs)
            elif type(a) is str:
                f.write(a)
            else:
                print('Couldn''t write ' + str(a))
        os.chdir(cwd)
        
    
    def reset(self, verbose=1):
        '''
        Retrives data for new object from lines read from file or resets 
        attribute values to those originally in the file
        '''
        print('loading attributes for this ' + self.name + ' molecule fresh from original file.')
        dictionary = lines_to_dictionary(self.file_lines)
        for (key, value) in dictionary.items():
                if 'calibration' in key:
                    self.add_calibration(value)
                else:
                    setattr(self, key, value)
         
    def write_new(self, f):
        for (attr, status) in self.attr_status.items():
            string = input('Enter ' + attr + ' for ' + self.name + ' as a float or anything else to skip.\n')
            if status == 0:             
                try: 
                    value = float(string)
                    self.attr_status = 2
                    setattr(self, attr, value) 
                    f.write(attr + '\t=\t' + str(value) + '\n')
                except ValueError:
                    print('skipped that for now.')
               # self.file_lines = f.readlines() #doesn't work. But what if I need to reset later?
          #  else:
          #      f.write(attr + '\t=\t' + str(self.attr) + '\n') #not necessary...
    
    def read_spectrum(self, lines):
        print('the read_spectrum function is not implemented yet')
        pass

    
    def add_calibration(self, calibration, useit=True, primary=True, writeit=False):
        '''
        'calibration' is a dictionary containing the calibration factor 'F_cal', 
        in C/mol; the mass measurement 'mass' for which it applies, as well as
        data about how it was obtained. This function adds the calibration to
        a list of calibration dictionaries for this instance of Molecule.
        If 'useit', 'F_cal' is used to set attribute 'F_<mass>' 
        If 'primary', this is the primary mass for measurement, and
        calibration['F_cal'] is (also) used to set 'F_cal' for this
        instance of Molecule for easier access, e.g. CO2.F_cal 
        If 'writeit', the calibration is written to the molecule's data file.
        '''
        self.calibrations += [calibration]   
        mass = calibration['mass']
        F_cal = calibration['F_cal']
        title = calibration['title']
        print('added calibration ' + title + ' for ' + self.name)
        if useit:       #use this calibration for this mass
            self.cal[mass] = calibration['F_cal']
            attribute_name = 'F_' + mass
            setattr(self, attribute_name, F_cal)
        if primary:     #use this mass by default
            self.primary = mass
            setattr(self, 'F_cal', F_cal)
        if writeit:     #write this calibration to file
            self.write(self.write_calibration, calibration)
    
    
    def new_calibration(self, calibration=None,
                        mass=None, F_cal=None, cal_type=None, 
                        chip=None, settings=None, notes=None,                    
                        expdate=None, andate=None, title=None, add_dates=True, 
                        useit=True, primary=True, writeit=True,
                        ):
        '''
        Puts values in a calibration dictionary and calls add_calibration
        '''
        andate = date_scott(andate)
        expdate = date_scott(expdate)        
        if title is None:
            title = expdate + '_' + andate
        elif add_dates:
            title = expdate + '_' + andate + '_' + title
        if chip == 'date':      
            #but ideally chip points to an object of the yet-to-be-written class chip
            chip = expdate
        if type(mass) is int:
            mass = 'M' + str(mass)
        if type(settings) is list:
            settings = {'SEM voltage':settings[0], 'speed':settings[1], 'range':settings[2]}
        if type(notes) is str:
            notes = [notes]
        calibration_i = {'title': title, 'F_cal': F_cal, 'mass': mass, 
                         'type': cal_type, 'experiment date': expdate, 
                         'analysis date': andate, 'chip': chip, 
                         'QMS settings': settings, 'Notes': notes}
        if calibration is None:
            calibration = calibration_i
        else:
            for (key, value) in calibration_i.items():
                if value is not None:
                    calibration[key] = value
                
        self.add_calibration(calibration, useit=useit, primary=primary, writeit=writeit)


    def read_calibration(self, lines, **kwargs):
        '''
        Generates a calibration dictionary from lines of text, which would 
        come from the molecule's data file. Then calls add_calibration.
        Never used anymore because reset does the same thing.
        '''
        calibration = lines_to_structure(lines)
        
        self.add_calibration(calibration, **kwargs)
        return calibration     
        

    def write_calibration(self, f, calibration): 
        '''
        Writes a calibration dictionary to the molecule's data file (pre-opened
        as 'f') in a way readable by read_calibration.
        Typically called by add_calibration.
        There's a much cleverer way to write this type of function.
        '''
        print('\nWriting calibration for ' + self.name + '\n')       
        
        title_line = 'calibration_' + calibration['title']
        lines = structure_to_lines(calibration, preamble=title_line)
        
        f.writelines(lines)      
        print('wrote calibration ' + calibration['title'] + ' for ' + self.name)
            
   
    def calibration_fit(self, mass='primary', ax1='new', useit=True, primary=True):
        if mass == 'primary':
            mass = self.primary
        title = mass + ' calibrations for ' + self.name
        n_mol = [0,]        #include 0,0 as a calibration point!
        Q_QMS = [0,]        #include 0,0 as a calibration point!
        for calibration in self.calibrations:
            if calibration['mass'] == mass:
                n_mol += [calibration['n_mol']]
                Q_QMS += [calibration['Q_QMS']]
        n_mol = np.array(n_mol)
        Q_QMS = np.array(Q_QMS)
        N = len(n_mol)
        pf1 = np.polyfit(n_mol, Q_QMS, 1)
        print('y = ' + str(pf1[0]) + ' x + ' + str(pf1[1]))
        F_cal = pf1[0]
        print()
        if useit:       #use this calibration for this mass
            attribute_name = 'F_' + mass
            print('using a fit value for ' + attribute_name + ' based on ' + str(N) + ' experiments.')
            self.cal[mass] = F_cal
            setattr(self, attribute_name, F_cal)
            if primary:
                self.F_cal = F_cal
                self.primary = mass
            
        pf_fun = np.poly1d(pf1)
        pf_x = np.array([min(n_mol), max(n_mol)])        
        
        if ax1 == 'new':
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
        if ax1 is not None:

            ax1.set_title(title)
            ax1.plot(n_mol*1e9, Q_QMS*1e9, 'k.', markersize=15)
            ax1.plot(pf_x*1e9, pf_fun(pf_x)*1e9, 'r--')
            ax1.set_xlabel('amount produced / nmol')
            ax1.set_ylabel('int. signal / nC')
   
    def set_temperature(self, T):
        self.T = T
        print('The set_temperature function is not implemented yet.')
        pass





if __name__ == '__main__':
    os.chdir('/home/soren/Dropbox (Soren Scott Inc)/Sniffer_Experiments/03_Pt_Sputtered/Analysis/16J14_time_response')
     
    if os.path.split(os.getcwd())[1] == 'EC_MS':   
                                #then we're running from inside the package
        data_directory = os.getcwd() + os.sep + 'data' + os.sep
    else:
        data_directory = os.getcwd() + os.sep + 'EC_MS' + os.sep + 'data' + os.sep    
    
    mol = Molecule('H2')    
    