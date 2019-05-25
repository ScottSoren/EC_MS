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
from numbers import Number

if os.path.split(os.getcwd())[1] == 'EC_MS':
                                #then we're running from inside the package
    import Chem
    from Object_Files import structure_to_lines, lines_to_dictionary
    from Object_Files import lines_to_structure, date_scott, update_lines
else:
    from . import Chem
    from .Object_Files import structure_to_lines, lines_to_dictionary
    from .Object_Files import lines_to_structure, date_scott, update_lines


preferencedir = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'preferences'
with open(preferencedir + os.sep + 'standard_colors.txt','r') as f:
    lines = f.readlines()
    standard_colors = lines_to_dictionary(lines, removecomments=False)['standard colors']
data_directory = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'data' + os.sep + 'molecules'
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
            self.attr_status['M'] = 0
            self.formula = name #but just put the name for now.
            self.M = Chem.get_mass(self.formula)
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
                    print('The file for ' + name + ' is empty!')
            self.reset(verbose=verbose)
            self.file_lines = ['name: ' + self.name] + self.file_lines
        except FileNotFoundError: # I don't know the name of the error, so I'll have to actually try it first.
            os.chdir(cwd)
            raise
            print('no file found for ' + name)
            if writenew:
                string = input('Start new data file for \'' + self.name + '\'? (y/n)\n')
                if string == 'n':
                    return None
                print('Writing new.\n')
                self.write(self.write_new)
            else:
                string = input('Start Molecule object \'' + self.name + '\' from scratch? (y/n)')
                if string == 'n':
                    return None
        os.chdir(cwd)

        if self.M == 0:
            print('setting self.M from formula') # debugging
            self.M = Chem.get_mass(self.formula)
        if not hasattr(self, 'molecule_mass'):
            print('setting self.molecule_mass from self.M') # debugging
            self.molecule_mass = self.M/Chem.NA

        self.transmission_function = None
        #self.color = self.get_color()


    def write(self, a=None, attr=None, *args, **kwargs):
        '''

        ... could move this to Object_Files, but a bit tricky

        this is supposed to be a versitile tool for writing to the Molecule's
        data file. Whether the added intricacy will be worth the lines of code
        it saves, I can't tell yet.

        Writes in one of the following ways:

        1. If the name of an attribute is given, that attribute is written.
        2. If a is a string, simply write a to the molecule's datafile.
        3. If a is a function, then it's a function that writes what you want
           to the data file given in the first argument
        4. If a is a dictionary, list, or tuple, convert it to lines according to
           the system encoded in Object_Files.
        5. If a is not given but keyword arguments are, write **kwargs to the
           data file according to the system encoded in Object_Files.
        '''


        if attr is not None:
            a = (attr, getattr(self, attr))
        elif a is None:
            if len(kwargs) == 0:
                print('nothing to write.')
                return
            else:
                a = kwargs

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

    def p_vap(self, T=None):
        if T is None:
            if 'T' in dir(self) and self.T is not None:
                T = self.T
            else:
                T = 298.15
        elif 'T' not in dir(self) or self.T is None:
            self.T = T
        return Chem.p_vap(self.name, T)


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
            if not self.primary is None and len(self.calibrations) > 0:
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
                string = input('Enter ' + attr + ' for ' + self.name + ' or whitespace to get default.\n')
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

    def get_RSF(self, RSF_source='NIST', mass='primary',
                transmission_function=None, verbose=True):
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
        if transmission_function is None:
            if self.transmission_function is None:
                print('WARNING: no transmission function for ' + self.name + '. using constant')
                def transmission_function(M):
                    return 1
            else:
                transmission_function = self.transmission_function #T(M)=1 unless otherwise stated
        elif self.transmission_function is None:
            self.transmission_function = transmission_function


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
                if verbose:
                    print(str(M) + ' ' + str(value))
                if M == 'Title':
                    continue
                self.RSF[M] = value / self.spectrum[mH] * vH
            if 'primary' in dir(self):
                self.rsf = self.RSF[self.primary]
            if verbose:
                print('RSFs from Hiden rsf, adjusted to mass of measurement according' +
                      'to NIST spectrum for ' + self.name)


        elif RSF_source == 'NIST':
            self.RSF = dict((key, value*transmission_function(int(key[1:]))) for key, value in self.IFCS.items())

            N2_M28_RSF = 2.283 * transmission_function(28)
            #print('N2_M28_RSF = ' + str(N2_M28_RSF)) # debugging
            self.RSF = dict([(key, value/N2_M28_RSF) for (key, value) in self.RSF.items()])

            if 'primary' in dir(self):
                self.rsf = self.ifcs * transmission_function(int(self.primary[1:])) #/ N2_M28_RSF
            if verbose:
                print('RSFs from ionization-fragmentation cross section in Ang^2 \'ifcs\'' +
                      'based on NIST cross section and spectrum for ' + self.name +
                      ' and the given transmission function T(M)')

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
        ax.bar(x, y)
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

    def mass_transfer_coefficient(self, system='chip', T=298.15,
                                  phi=0.5, l_p=100e-6, d_p=20e-9,
                                  p_chip=1e5, n_dot_0=2.6e-9, A=0.196e-4,
                                  verbose=True):
        K_H = self.kH * Chem.R * T
        self.K_H = K_H
        if verbose:
            print('K_H = ' + str(K_H*1e-2) + ' bar/mol')
        if system == 'chip':
            self.h = K_H * n_dot_0 / (p_chip * A)
        else:
            self.h = K_H * phi * d_p / (3 * l_p) \
                     * np.sqrt(8 / (np.pi * Chem.R * T * self.M*1e-3))
        if verbose:
            print('h = ' + str(self.h) + ' m/s')
        return self.h

    def set_temperature(self, T):
        self.T = T
        print('The set_temperature function is not implemented yet.')
        pass

    def get_bg(self, *args, **kwargs):
        x, y = self.get_flux(*args, **kwargs, removebackground=False)
        background = np.mean(y)
        self.background = background
        return background

    def get_color(self):
        try:
            return self.color
        except AttributeError:
            try:
                return standard_colors[self.primary]
            except AttributeError:
                print('WARNING: ' + str(self) + ' has no attribute \'primary\'')
            except KeyError:
                print('WARNING: standard_colors has no entry for ' + str(self.primary))

    def get_flux(self, MS_data, tspan='tspan', density=None,
                 unit='pmol/s', verbose=True, override=False, x=None,
                 removebackground=None, background=None, t_bg=None,
                 endpoints=5):
        '''
        returns [x, y] where x is the t corresponding to the primary mass of the
        molecule in 'mol' and y is the molecular flux in nmol/s, calculated from
        the MS_data for its primary mass and the value of F_cal read from the
        molecule's text file.
        '''
        if verbose:
            print('calculating flux of ' + self.name)

        if hasattr(self, 'cal_mat'):
            cal_mat = self.cal_mat
        else:
            print('no cal_mat! using self.primary and self.F_cal instead.')
            F_cal = self.F_cal
            mass = self.primary
            cal_mat = {mass:1/F_cal}

        if tspan is None:
            if x is not None:
                tspan = [x[0], x[-1]]
            else:
                tspan = 'tspan'
        if type(tspan) is str and not tspan=='all':
            tspan = MS_data[tspan]

        if x is None:
            if density is None:
                x = MS_data[self.primary + '-x']
                mask = np.logical_and(tspan[0]<x, x<tspan[-1])
                # Don't cut off outer endpoints before evt interpolation (if used by plot_vs_potential)
                extra_left = np.append(mask[1:], False)
                extra_right = np.append(False, mask[:-1])
                mask = np.logical_or(extra_left, extra_right)
                x = x[mask]
            else:
                x = np.linspace(tspan[0], tspan[-1], density*np.floor(tspan[-1] - tspan[0]))
            # ^ approx 5 datapoints a second

        y = 0
        for mass, C in cal_mat.items():
            x0 = MS_data[mass + '-x']
            s0 = MS_data[mass + '-y']
            s = np.interp(x, x0, s0)
            y += s * C # units: [A] * [mol/C] = [mol/s]

        if (background is None or background is 'preset') and hasattr(self, 'background'):
            background = self.background
        elif background == 'preset':
            background = None

        if removebackground is None and background is not None:
            removebackground = True

        if removebackground:
            if background is None:
                background = 'constant'
            if background in ['start', 'begining', 'first']:
                background = np.mean(y[:endpoints])
            elif background in ['finish', 'end', 'last']:
                background = np.mean(y[-endpoints:])
            elif background=='constant':
                print('using minimum value as constant background for ' + self.name)
                if type(removebackground) is float:
                    background = removebackground * min(y)
                elif t_bg is not None:
                    try:
                        mask = np.logical_and(t_bg[0]<x, x<t_bg[-1])
                        background = np.mean(y[mask])
                    except TypeError:
                        background = np.interp(t_bg, x, y)
                else:
                    background = min(y)
                if not hasattr(self, 'background') or self.background is None:
                    self.background = background
            elif background=='linear':
                x_end = [np.average(x[:endpoints]), np.average(x[-endpoints:])]
                y_end = [np.average(y[:endpoints]), np.average(y[-endpoints:])]
                background = np.interp(x, x_end, y_end)
                print('using linear background for ' + self.name)
            elif isinstance(background, Number):
                #background = background
                if verbose:
                    print('using preset constant background for ' + self.name)
                pass
            y = y - 0.99*background #so that we don't break the log scale.
            #I should get rid of this and assume the caller knows what they're doing.

        # important that this comes after the background bit, so that we don't
        # subtract a background in different units than the signal
        if 'nmol' in unit:
            y = y * 1e9
        elif 'pmol' in unit:
            y = y * 1e12
        if 'cm$^{-2}$' in unit or '/cm^2' in unit:
            y = y / MS_data['A_el']

        return x, y





def reset_datafiles(mols, attrs, mdict={}):
    '''
    loads all of the molecules in mols and rewrites their data files with
    only the attributes listed in attrs.
    Look through the .git history if you need any old information removed by
    this function.
    Tuple entries attrs[i] rename the attribut from attrs[i][0] to attrs[i][1]
    '''
    for mol in mols:
        if mol in mdict:
            continue
        elif type(mol) is str:
            mdict[mol] = Molecule(mol) #might return None.
        else: #then mol is a Molecule object and mol.name is the key
            mdict[mol.name] = mol

    cwd = os.getcwd()
    os.chdir(data_directory)
    for mol, m in mdict.items():
        if m is None:
            continue

        f = open(mol +'.txt', 'w')
        f.write('#cleaned up ' + date_scott() + '\n')
        f.close()
        print('wiped ' + mol + '.txt clean.')
        #I think this makes blank datafiles that the molecules can then
        #write to.
        for attr in attrs:
            if type(attr) is list or type(attr) is tuple: #rename attr
                var = attr[0]
                newvar = attr[1]
            else:
                var = attr
                newvar = attr
            try:
                value = getattr(m, var) #keep the attribute's name
                m.write((newvar, value))
            except AttributeError:
                print('no attribute ' + var + ' for Molecule ' + mol)

    os.chdir(cwd)
    return mdict



def add_to_datafiles(attr, d, mdict={}, mols='all'):
    '''
    d is a dictionary containing attribute attr for specified molecules.
    This function writes that attribute to each datafile
    '''
    for (key, value) in d.items():
        if type(key) is str:
            if not mols=='all' and key not in mols:
                continue
            if key in mdict:
                m = mdict[key]
            else:
                m = Molecule(key)
                mdict[key] = m
        if m is not None:
            setattr(m, attr, value)
            m.write((attr, value))
    return mdict



def add_script_to_datafiles(path, file_name, attrs='all', mdict={}, mols='all'):
    '''
    sorts data stored in another script in dictionary form into molecule data
    form.
    Tuple entries attrs[i] rename the attribut from attrs[i][0] to attrs[i][1]
    '''
    module_name = file_name.split('.')[0] #drop extension
    cwd = os.getcwd()
    os.chdir(path)
    module = __import__(module_name)
    os.chdir(cwd)

    check = False
    if attrs == 'all':
        check = True
        attrs = (d for d in dir(module) if type(d) is dict)

    for attr in attrs:
        if check: #then check manually
            var = attr
            newvar = input('Write data from \'' + var + '\'? (y/<newname>/n)\n')
            if newvar == 'n':
                continue
            elif newvar == 'y':
                newvar = var
        elif type(attr) is list or type(attr) is tuple:
            var = attr[0]
            newvar = attr[1]
        else:
            var = attr
            newvar = attr

        d = getattr(module, var)
        mdict = add_to_datafiles(newvar, d, mdict, mols=mols)




if __name__ == '__main__':


    plt.close('all')


    reset17I21 = True
    if reset17I21:
        #clean up old useless calibrations etc,
        #add new data to data files, from Daniel's script

        mols = ['Ar', 'C2H4', 'C2H6', 'CH4', 'N2','O2',
                'Cl2', 'CO', 'CO2', 'H2', 'H2O', 'He',
                'acetaldehyde', 'ethanol',
                'air','Henriksen2009',
                ]

        mdict = dict((mol, Molecule(mol)) for mol in mols)


        #CLEANUP:

        keep = ['name', 'formula', 'M', 'D', 'kH', #'primary',
                'sigma_100eV', 'spectrum', 'RSF_Hiden', ('Hiden', 'RSF_Hiden')]

        reset_datafiles(mols=mols, attrs=keep, mdict=mdict)


        #ADD NEW DATA:
        file_name = 'Membrane_chip.py'
        path = '/home/scott/Dropbox/python_master/data'

        new = (('m', 'molecule_mass'), ('s', 'molecule_diameter'),
               ('eta', 'dynamic_viscosity'), ('rho', 'density_RTP'),
               ('D', 'D_gas_RTP'))

        add_script_to_datafiles(path, file_name, attrs=new, mdict=mdict)




    testtxt = False
    if testtxt:
        H2 = Molecule('H2')
        print(H2.F_cal)
        H2.F_cal = 1.111
        H2.rewrite('data/test.txt')
        test = Molecule('test')
        print(test.F_cal)
        test.calibration_fit()


    plotcals=False
    if plotcals:
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


    if True: # write primary to molecules based on spectrum
        mols = ['Ar', 'C2H4', 'C2H6', 'CH4', 'N2','O2',
                'Cl2', 'CO', 'CO2', 'H2', 'H2O', 'He',
                'acetaldehyde', 'ethanol',
                'air','Henriksen2009',
                ]

        for mol in mols:
            m = Molecule(mol)
            try:
                mass, value = zip(*[(k, v) for (k, v) in m.spectrum.items() if type(v) in [int, float]])
            except AttributeError:
                continue
            if len(value)>0:
                i = np.argsort(np.array(value))
                primary = mass[i[-1]]
                print(mol + ' ' + primary)

                m.primary = primary
                m.write(attr='primary')

