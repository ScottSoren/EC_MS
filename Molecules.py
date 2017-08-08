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
    from Object_Files import structure_to_lines, lines_to_dictionary 
    from Object_Files import lines_to_structure, date_scott, update_lines
else:
    from . import Chem
    from .Object_Files import structure_to_lines, lines_to_dictionary 
    from .Object_Files import lines_to_structure, date_scott, update_lines

data_directory = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'data'
cwd = os.getcwd()
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
    def __init__(self, name, formula=None, writenew=True, verbose=True):
        self.name = name
        self.cal = {}
        self.__str__ = '<' + name + ', instance of EC_MS class \'Molecule\'>'
        self.primary = None     #the primary mass to measure at
        self.calibrations = []      #will store calibration data
        self.attr_status = {'D':0,'kH':0,'n_el':0}
        if formula is None:
            self.attr_status['formula'] = 0 #so that it asks me for a formula when initializing the file.
            self.formula = name #but just put the name for now.
        else:
            self.formula = formula
        # 0 for undefined, 1 for loaded from file, 2 for set by function
        file_name = self.name + '.txt'
        cwd = os.getcwd()
        os.chdir(data_directory)
        try: 
            with open(file_name,'r') as f:
                self.file_lines = f.readlines()
                if len(self.file_lines) == 0:
                    print('The file for ' + name + ' is empty! Writing new.')
                    raise FileNotFoundError
            self.reset(verbose=verbose)
            self.file_lines = ['name: ' + self.name] + self.file_lines
        except FileNotFoundError: # I don't know the name of the error, so I'll have to actually try it first.
            print('no file found for ' + name)  
            if writenew:
                print('Writing new.\n')
                self.write(self.write_new)
        os.chdir(cwd)
        
        self.M = Chem.Mass(self.formula) #molar mass from chemical formula
    
    
    def write(self, a, *args, **kwargs):
        cwd = os.getcwd()
        os.chdir(data_directory)
        file_name = self.name + '.txt'
        with open(file_name, 'a') as f:
            if callable(a):
                a(f, *args, **kwargs)
            elif type(a) is str:
                if a[-1] != '\n':
                    a += '\n'
                f.write(a)
            elif type(a) in (list, dict):
                if 'key' in kwargs.keys():
                    lines = structure_to_lines(a, preamble=kwargs['key'])
                else:
                    lines = structure_to_lines(a)
                f.writelines(lines)
            elif type(a) is tuple:
                lines = structure_to_lines(a[1], preamble=a[0])
                f.writelines(lines)
            else:
                print('Couldn''t write ' + str(a))
        os.chdir(cwd)
    
    
    def rewrite(self, file='default'):
        if file=='default':
            file = self.name + '.txt'
        
        newlines = update_lines(self.file_lines, self.__dict__, 
                                oldkeys=['file_lines', 'calibrations', 'attr_status', '__str__'],
                                )
        
        if type(file) is str:
            os.chdir(data_directory)
            with open(file, 'w') as f:
                f.writelines(newlines)
            os.chdir(cwd)
        else: 
            file.writelines(newlines)
        
    
    def reset(self, verbose=True):
        '''
        Retrives data for new object from lines read from file or resets 
        attribute values to those originally in the file
        '''
        if verbose:
            print('loading attributes for this ' + self.name + ' molecule fresh from original file.')
        dictionary = lines_to_dictionary(self.file_lines)
        for (key, value) in dictionary.items():
            if 'calibration' in key:
                self.add_calibration(value)
            if 'Spectrum' in key:
                self.spectrum = value
            else:
                setattr(self, key, value)
        
        if 'F_cal' not in dir(self) and 'primary' in dir(self):
            if not self.primary is None:
                self.F_cal =  self.calibration_fit(mass='primary', ax=None,
                                                   useit=True, primary=True, verbose=True)
        elif type(self.F_cal) is list and len(self.F_cal) == 1:
            self.F_cal = self.F_cal[0]
            #print('F_cal was list of length 1. Now, F_cal = ' + str(self.F_cal))
        else:
            #print('F_cal was as it should be')
            pass
                
         
    def write_new(self, f):
        for (attr, status) in self.attr_status.items():
            if status == 0: 
                string = input('Enter ' + attr + ' for ' + self.name + ' or just whitespace.\n')
                if len(string.strip()) == 0:
                    print('skipped that for now.')
                    continue
                try: 
                    value = float(string)
                    self.attr_status = 2
                    setattr(self, attr, value) 
                    f.write(attr + '\t=\t' + str(value) + '\n')
                except ValueError:
                    print('not a float but okay.')
                    value = string
                    self.attr_status = '2' #just for irony, I'll save this status as not a float.
                    setattr(self, attr, value) 
                    f.write(attr + ': ' + str(value) + '\n')                    
                    

               # self.file_lines = f.readlines() #doesn't work. But what if I need to reset later?
          #  else:
          #      f.write(attr + '\t=\t' + str(self.attr) + '\n') #not necessary...
    
    def get_RSF(self, RSF_source='NIST', mass='primary'):
        '''
        Requires that a spectrum and total ionization cross section are already 
        loaded, and preferably also a relative sensitivity factor.
        Generates dictionaries of ionization-fragmentation cross sections 'IFCS' 
        and 'RSF' for each mass in the spectrum. Saves the respective value for
        the stated primary mass also as 'ifcs' and 'rsf'.
        
        '''
        self.IFCS = {}       #this will be a dictionary containing the electron ionization
                            #portion of the 'relative sensitivity' for each mass,
                            #i.e. the partial ionization cross-section in Ang^2 at 100keV
        spec_total = 0
        for value in self.spectrum.values():
            if type(value) is not str:
                spec_total += value
                
        for (M, value) in self.spectrum.items():
            if M == 'Title':
                continue
            self.IFCS[M] = value / spec_total * self.sigma_100eV
        if 'primary' in dir(self):
            self.ifcs = self.IFCS[self.primary]
            
        if RSF_source == 'Hiden':
            try:
                self.RSF = {}       #this will be the relative sensivity factor at each
                        #mass, where N2 at M28 is 1, from Hiden Analytical  
                mH = self.Hiden[0]
                vH = self.Hiden[1]
            except AttributeError: 
                print('no Hiden RSF found for ' + self.name)
                return None
            for (M, value) in self.spectrum.items():
                print(str(M) + ' ' + str(value))
                if M == 'Title':
                    continue
                self.RSF[M] = value / self.spectrum[mH] * vH
            if 'primary' in dir(self):
                self.rsf = self.RSF[self.primary]
            print('returning Hiden rsf, adjusted to mass of measurement according' +
                  'to NIST spectrum for ' + self.name)                    
                    
                    
        elif RSF_source == 'NIST':
            self.RSF = self.IFCS
            if 'primary' in dir(self):
                self.rsf = self.ifcs
            print('returning ionization-fragmentation cross section in Ang^2 \'ifcs\'' +
                  'based on NIST cross section and spectrum for ' + self.name)
            
        if mass == 'primary':
            mass = self.primary
        
        return self.RSF[mass]

    def plot_spectrum(self, top=100, ax='new'):
        if ax == 'new':
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)
        x = []
        y = []
        for (mass, value) in self.spectrum.items():
            if mass == 'Title':
                continue
            x += [int(mass[1:])]
            y += [value]
        y = np.array(y) / max(y) * top
        x = np.array(x)
        ax.bar(x - 1/2, y)
        ax.set_xticks(x)
        ax.set_xticklabels([str(m) for m in x])
        ax.set_title('literature QMS spectrum for ' + self.name)
        
    
    def add_calibration(self, calibration, 
                        useit=True, primary=True, writeit=False, verbose=False):
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
        if verbose:
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
            
   
    def calibration_fit(self, mass='primary', ax='new', color='k', plotfactor=1,
                        useit=True, primary=True, verbose=True):
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
        if verbose:
            print('y = ' + str(pf1[0]) + ' x + ' + str(pf1[1]))
        F_cal = pf1[0]
        print()
        if useit:       #use this calibration for this mass
            attribute_name = 'F_' + mass
            if verbose:
                print('using a fit value for ' + attribute_name + ' based on ' + str(N) + ' experiments.')
            self.cal[mass] = F_cal
            setattr(self, attribute_name, F_cal)
            if primary:
                self.F_cal = F_cal
                self.primary = mass
            
        pf_fun = np.poly1d(pf1)
        pf_x = np.array([min(n_mol), max(n_mol)])        
        
        if ax == 'new':
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)
        if ax is not None:

            ax.set_title(title)
            ax.plot(n_mol*1e9*plotfactor, Q_QMS*1e9*plotfactor, '.', color=color, markersize=15)
            ax.plot(pf_x*1e9*plotfactor, pf_fun(pf_x)*1e9*plotfactor, '--', color=color,)
            ax.set_xlabel('amount produced / nmol')
            ax.set_ylabel('int. signal / nC')
        return F_cal
   
    def set_temperature(self, T):
        self.T = T
        print('The set_temperature function is not implemented yet.')
        pass





if __name__ == '__main__':
    plt.close('all')
    
    H2 = Molecule('H2')
    print(H2.F_cal)
    H2.F_cal = 1.111
    H2.rewrite('data/test.txt')
    test = Molecule('test')
    print(test.F_cal)
    test.calibration_fit()
    
    
    
    mols = {'CO2':('brown', 10),'H2':('b', 1),'O2':('k', 1)}
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for (mol, (color, plotfactor)) in mols.items():
        m = Molecule(mol)
        m.calibration_fit(ax=ax, color=color, plotfactor=plotfactor)
    ax.set_title('')
    #plt.savefig('Internal_calibrations.png')
    
    

    #mol = Molecule('C2H4')  
    #mol.plot_spectrum()
    #mol = Molecule('C2H6')  
    #mol.plot_spectrum()
    
    
    
    
    
    
    