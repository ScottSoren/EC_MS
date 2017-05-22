# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 13:44:31 2016

@author: scott
"""

from __future__ import print_function, division
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

if os.path.split(os.getcwd())[1] == 'EC_MS':   
                                #then we're running from inside the package
    from Molecules import Molecule
    from Combining import cut
    from Object_Files import lines_to_dictionary, date_scott
    from EC import sync_metadata
    import Chem
    data_directory = os.getcwd() + os.sep + 'data' + os.sep
else:
    from .Molecules import Molecule
    from .Combining import cut
    from .Object_Files import lines_to_dictionary, date_scott
    from .EC import sync_metadata
    from . import Chem
    data_directory = os.getcwd() + os.sep + 'EC_MS' + os.sep + 'data' + os.sep



def rewrite_spectra(NIST_file='default', RSF_file='default',
                    writesigma=False, writespectra=False, writeRSF=False):
    '''
    Reformats NIST data copied to a text file given by 'file', and writes it to 
    the molecule files in data_directory
    '''
    if writesigma or writespectra:
        if NIST_file == 'default':
            NIST_file = data_directory + 'NIST_spectra_data.txt'
        with open(NIST_file, 'r') as f:
            lines = f.readlines()
        structure = lines_to_dictionary(lines)
        sigma_dict = structure['sigma_100eV'] # dictionary with electron impact ionizations in Ang^2
        spectra_dict = structure['Spectra'] # dictionary with QMS spectra as mass1,rel_val1 mass2,relval2 
    if writesigma:
        for (mol, sigma) in sigma_dict.items():
            m = Molecule(mol)
            l = ('sigma_100eV', sigma)
            m.write(l)
    if writespectra:
        for (mol, specline) in spectra_dict.items():
            m = Molecule(mol)
            masses = []
            rel_vals = []        
            spectrum = {'Title': 'NIST_16L14'}
            for spec in specline.split(' '):
                (mz, v) = spec.split(',')
                masses += ['M' + mz]
                rel_vals += [eval(v)]
            rel_vals = np.array(rel_vals)
            rel_vals = 100 * rel_vals / np.max(rel_vals)       #normalize spectrum to % of max peak
            for (mass, rel_val) in zip(masses, rel_vals):
                spectrum[mass] = rel_val               
            l = ('Spectrum', spectrum)
            m.write(l)      
    if writeRSF:
        if RSF_file == 'default':
            RSF_file = data_directory + 'Relative_Sensativity_Factors.txt'
        with open(RSF_file) as f:
            lines = f.readlines()
        structure = lines_to_dictionary(lines)
        for (mol, value) in structure.items():
            m = Molecule(mol)
            mass = 'M' + str(value[0])
            RSF = value[1]
            l = ('Hiden', (mass, RSF))
            m.write(l)
    return structure

    
def calibration_compare(calmol = ['H2','O2','CO2','Cl2'] ):
    '''
    Checks the calibration factors calculated internally at the sniffer setup
    against predicted sensitivities based on the NIST database. 
    There is no correlation. We are crazy sensitive to H2.
    '''
    for mol in calmol:
        m = Molecule(mol, verbose=False)
        m.get_RSF()
        m.calibration_fit(verbose=False, ax=None)
        p = m.primary
        c1 = m.F_cal
        c2 = m.ifcs

        r12 = c1/c2
        RSFit = 'rsf' in dir(m)
        if RSFit:
            c3 = m.rsf
            r13 = c1/c3
            
        print(mol + ' at ' + p + ', calibration factors:\n' +
            '\tc1 = ' + str(c1) + ' C/mol (experimental from Pt) \n' +
            '\tc2 = ' + str(c2) + ' Ang^2 (calculated from NIST) \n' + 
            'r12 = ' + str(r12) + '  (ratio) ')
        if RSFit:
            print('\tc3 = ' + str(c3) + '   (relative sensitivity factor from Hiden) \n' + 
                'r13 = ' + str(r13) + '  (ratio)')          
            
            
def RSF_to_F_cal(quantmol = {'H2':'M2', 'He':'M4', 'CH4':'M15', 'H2O':'M18',
                             'N2':'M28', 'CO':'M28', 'C2H4':'M26', 'C2H6':'M30',
                             'O2':'M32', 'Ar':'M40', 'CO2':'M44', 'Cl2':'M70'},
                 calmol = {'O2':'M32', 'CO2':'M44'},
                 ax = 'new', writeit = False, writeprimary=False):
    '''
    Generates calibration factors 'F_cal' for primary masses for all molecules
    in 'quantmol' that don't already have one. The new F_cal are based on the 
    internal calibrations for a given set of molecules given by 'calmol' and the
    Hiden RSF values. F_cal is assumed proportional to RSF, with r = F_cal/RSF.
    
    'primary' and 'F_cal' are written to each molecule's file.
    
    ----
    As of 16L15, only CO2 and O2 should be used to calculate r, as we are
    freakishly sensitive to H2, there is no Hiden RSF for Cl2, and we have no
    other reliable internal calibrations. We need more internal calibrations.
    '''
    RSFs = []
    F_cals = []
    RSF_dict = {}
    for (mol, mass) in calmol.items():
        m = Molecule(mol, verbose=False)
        F_cal = m.calibration_fit(mass=mass, ax=None, useit=True, primary=True)
        rsf = m.get_RSF()
        F_cals += [F_cal]
        RSFs += [rsf]
    
    def fit_fun(x, a):
        return a*x
    r, pcov = curve_fit(fit_fun, RSFs, F_cals, p0 = 10)
    
    print('r = ' + str(r) + 'C/mol')
    
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    for (mol, mass) in quantmol.items():
        m = Molecule(mol, verbose=False)
        m.primary = mass
        rsf = m.get_RSF()
        if rsf is None:
            print('missing rsf for ' + mol)
            continue
        RSF_dict[mol] = rsf
        if len(m.cal) > 0:      #prioritize our internal calibrations!
            F_cal = m.calibration_fit(mass=mass, ax=None, useit=True, primary=True)
            ax.plot(rsf, F_cal, '.r', markersize=15)
            if writeit:
                m.write('#fit to calibrations on mass ' + mass + ' ' + date_scott())
                l = ('F_cal', F_cal)
                m.write(l)
            continue
        F_cal = r * rsf
        ax.plot(rsf, F_cal, '.b', markersize=15)
        if writeit:
            m.write('#extrapolated from rsf on mass ' + mass + ' ' + date_scott())
            l = ('F_cal', F_cal)
            m.write(l)
        if writeprimary:
            l = ('primary', mass)
            m.write(l)
        print(mol + ': F_cal = ' + str(F_cal))
            
    ax.set_xlabel('Relative Sensitivity Factor')
    ax.set_ylabel('F_cal / [C/mol]')
    tickmol = ['C2H6', 'C2H4', 'CH4', 'H2', 'CO2']
    ticknr = [0, 1, 1.5]
    #ax.set_xticks(ticknr + [RSF_dict[i] for i in tickmol])
    #ax.set_xticklabels([str(n) for n in ticknr] +tickmol)
    ax.set_ylim([0,20])


def get_flux(MS_data, mol, tspan='tspan_2', removebackground=False, 
             unit='pmol/s', verbose=True):
    '''
    returns [x, y] where x is the t corresponding to the primary mass of the
    molecule in 'mol' and y is the molecular flux in nmol/s, calculated from 
    the MS_data for its primary mass and the value of F_cal read from the 
    molecule's text file.
    '''    
    if type(mol) is str:
        m = Molecule(mol, verbose=False)
    else:
        m = mol
    if verbose:
        print('calculating flux of ' + m.name)
    
    mass = m.primary
    F_cal = m.F_cal
    x = MS_data[mass + '-x']
    y = MS_data[mass + '-y'] / F_cal # units: [A]/[C/mol]=[mol/s]
    if 'nmol' in unit:
        y = y * 1e9
    elif 'pmol' in unit:
        y = y * 1e12
    if 'cm^2' in unit:
        y = y / MS_data['A_el']
    
    if type(tspan) is str:
        tspan = MS_data[tspan]
    if tspan is not None:
        x, y = cut(x,y,tspan) 
    
    if removebackground:
        if type(removebackground) is float:
            background = removebackground * min(y)
        else:
            background = min(y)
        y = y - 0.99*background
    return [x,y]    
            
def predict_current(EC_and_MS, mols, tspan=None, RE_vs_RHE=None, A_el=None,
                    ax='new', colors=None, verbose=1):
    '''
    calculates a predicted electrical current based on MS_data and molecular
    data loaded from the files named in the list 'mols.'
    As of 16K28 just uses the primary masses
    tspan = x1 uses the time axis of the first primary mass.
    '''        
    V_str, J_str = sync_metadata(EC_and_MS, RE_vs_RHE, A_el)
    A_el = EC_and_MS['A_el']
    if A_el is None:
        A_el = 1
    
    partials = []
    for mol in mols:
        molecule = Molecule(mol)
        '''
        mass = molecule.primary
        x = EC_and_MS[mass + '-x']
        y = EC_and_MS[mass + '-y']
        F_cal = molecule.F_cal
        i_mol = y / F_cal       #molecular flux in mol/s
        '''
        [x, i_mol] = get_flux(EC_and_MS, molecule)
        j_mol = i_mol / A_el    #molecular flux densith in mol/cm^2/s
        j = j_mol * molecule.n_el * Chem.Far  * 1e3      #current density in mA/cm^2
        partials += [[x,j]]
    
    #this handling of tspan might be excessive for now...
    #it's intricate because I might need to interpolate MS data to EC time
    if tspan == None and 'tspan_2' not in EC_and_MS.keys():
        tspan == 'x1:'
    if tspan == 'x1':
        t = partials[0][0]
    elif 'time/s' in EC_and_MS.keys():
        t = EC_and_MS['time/s']
    if tspan == None or tspan == 'tspan_2':
        tspan = EC_and_MS['tspan_2']
    if len(tspan) == 2 and type(tspan) is not str:
        t = [t_i for t_i in t if tspan[0]<t_i and t_i<tspan[1]]
    else:
        tspan = [t[0], t[-1]]
    
    if ax=='new':
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    js = []        
    j_total = np.zeros(np.shape(t))
    for ([x,j], mol) in zip(partials, mols):
        I_keep = [I for (I, x_I) in enumerate(x) if t[0]<x_I and x_I<t[-1]]
        x = x[I_keep]
        j = j[I_keep]
        j_int = np.interp(t, x, j)
        j_total += j_int
        js += [j_int]
        if ax is not None and colors is not None:
            ax.plot(t, j_int, colors[mol], label=mol)

    if ax is not None:
        ax.plot(t, j_total, 'k:', label='total')


         

if __name__ == '__main__':
    plt.close('all')
    #calibration_compare()
    RSF_to_F_cal(calmol = {'CO2':'M44'}, writeprimary=False)
    
    plt.savefig('Fcal_vs_RSF.png')
    
        
    
    
    
    