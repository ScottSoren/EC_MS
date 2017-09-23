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


preferencedir = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'preferences' 
with open(preferencedir + os.sep + 'standard_colors.txt','r') as f:
    lines = f.readlines()
    standard_colors = lines_to_dictionary(lines)['standard colors']


def rewrite_spectra(NIST_file='default', RSF_file='default', mols='all',
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
            if not (mols == 'all' or mol in mols):
                continue
            m = Molecule(mol)
            l = ('sigma_100eV', sigma)
            m.write(l)
    if writespectra:
        for (mol, specline) in spectra_dict.items():
            if not (mols == 'all' or mol in mols):
                continue
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
            if not (mols == 'all' or mol in mols):
                continue
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

            
def RSF_to_F_cal(*args, **kwargs):
    print('\'RSF_to_F_cal\' has been renamed \'recalibrate\'. Remember that next time!')
    return recalibrate(*args, **kwargs)


def recalibrate(quantmol = {'H2':'M2', 'He':'M4', 'CH4':'M15', 'H2O':'M18',
                             'N2':'M28', 'CO':'M28', 'C2H4':'M26', 'C2H6':'M30',
                             'O2':'M32', 'Ar':'M40', 'CO2':'M44', 'Cl2':'M70'},
                                #molecules we want to calc F_cal (at the given mass) for by extrapolation
                mdict = {},
                calibrate = {'CO2':'M44',},
                              #molecules to base the extrapolation on
                F_cals = {}, #new calibrations to base the extrapolation on. 
                             #if 
                RSF_source = 'NIST', #'NIST' is partial ionization cross section, citable!
                                      #'Hiden' is RSF relative to N2, uncitable!
                internal = {},
                external = {},
                transmission_function='default',
                trust = 'new',  #defines the red points, i.e.
                #which calibrations to trust and use over extrapolated values,
                #  trust='default': do a calibration fit and trust it for molecules with saved calibration entries
                #  trust='new': trust the F_vals explicitly given in 
                #  can also put a list.
                trendline = True,
                ax = 'new', 
                writeit = False, writeprimary=False, rewriteit=False):
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
    
    ----
    largely rewritten 17H08
    '''
    print('\n\nfunction \'recalibrate\' at your service!\n')


    RSF_unit = {'Hiden':'a.u.', 'NIST':'a.u.'}[RSF_source]

    molset = set(quantmol.keys()) | set(internal.keys()) | set(external.keys()) | set(F_cals.keys())   
    for mol in molset:
        if mol not in mdict:
            mdict[mol] = Molecule(mol)
    #mdict = dict((mol, Molecule(mol, verbose=False)) for mol in molset)
   
    if type(calibrate) is dict:
        calmol = calibrate
    else:
        calmol = {}
    if  calibrate == 'internal' or calibrate == 'all':
        calmol.update(dict((mol, mdict[mol].primary) for mol in set(internal.keys())))
    if calibrate == 'external' or calibrate == 'all':
        calmol.update(dict((mol, mdict[mol].primary) for mol in set(external.keys())))
    
    print(calmol)
    #trust = [mdict[mol] for mol in trust] #trust should be molecule objects, not strings. #FALSE
    
    F_cals.update(internal)
    F_cals.update(external)
    print('F_cals = ' + str(F_cals))

    if trust == 'default':
        trust = {mol for mol, m in mdict.items() if 'calibrations' in m.__dict__}
    elif trust == 'new':
        trust = set(F_cals.keys())

    
    if len(internal) == 0 and len(external) == 0:
        internal = calmol
    
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        
    RSF_vec = []
    F_cal_vec = []
    
    if transmission_function == 'default':
        def T(M):
            return np.pow(M, -1/2)
    elif transmission_function == 1:
        def T(M):
            return 1
    else:
        T = transmission_function


    for (mol, mass) in calmol.items() :
        F_cal = F_cals[mol]
        m = mdict[mol]
        #mass = m.primary
        m.F_cal = F_cal
        print('F_' + mol + '_' + mass + ' = ' + str(F_cal))
        rsf = m.get_RSF(RSF_source=RSF_source, transmission_function=T, mass=mass)
        F_cal_vec += [F_cal]
        RSF_vec += [rsf]
    
    def fit_fun(x, a):
        return a*x
  
    r, pcov = curve_fit(fit_fun, RSF_vec, F_cal_vec, p0 = 10)  
    print('r = ' + str(r) + ' (C/mol)/' + RSF_unit)
        
    RSF_dict = {}   
    for (mol, m) in mdict.items():
        if mol in quantmol:
            mass = quantmol[mol]
        elif mol in calmol:
            mass = calmol[mol]
        else:
            print('if you include an F_cal in F_cals, make sure to also ' +
                  'include the molecule in quantmols, or it\'s ignored.')
            mass = m.primary
            
        rsf = m.get_RSF(RSF_source=RSF_source, transmission_function=T, mass=mass)
        if writeit:
            m.write('#the folowing rsf is calculated for ' + mass + ' on ' + date_scott())
            m.write(('rsf', rsf))
        if rsf is None:
            print('missing rsf for ' + mol)
            continue
        RSF_dict[mol] = rsf
        if mol in trust:      #The red points.
            if mol in F_cals:
                F_cal = F_cals[mol]
            elif 'calibrations' in m.__dict__: #changed mol to m, 17I08
                try:
                    F_cal = m.calibration_fit(mass=mass, ax=None, useit=True, primary=True)
                except ValueError:
                    F_cal = m.F_cal
            elif 'cal' in m.__dict__:
                F_cal = m.cal[mass]
            elif 'F_cal' in m.__dict__:
                F_cal = m.F_cal
            else:
                print('No value of F_cal to trust for ' + mol + '!')  
                break #becomes blue.
            if mol in internal:
                ax.plot(rsf, F_cal, 's', color=standard_colors[mass], markersize=10)
            else: #F_cals not given in internal are plotted as if they were external calibrations
                ax.plot(rsf, F_cal, '^', color=standard_colors[mass], markersize=10)
            if writeit:                
                m.write('#the following F_cal value is for ' + mass + ', trusted on ' + date_scott())
                l = ('F_cal', F_cal)
                m.write(l)
        else:
            F_cal = r * rsf  #the extrapolation!
            ax.plot(rsf, F_cal, '.', color=standard_colors[mass], markersize=10)
        if 'ca' in m.__dict__:
            m.cal[mass] = F_cal
        if 'primary' not in m.__dict__:
            m.primary = mass
            m.F_cal = F_cal
        elif m.primary == mass:
            m.F_cal = F_cal
        
        if writeit:
            m.write('#the following F_cal value is for ' + mass + ', extrapolated ' +
                    'from trusted values based on RSF from ' + RSF_source + ' on ' + date_scott())
            l = ('F_cal', F_cal)
            m.write(l)
        if writeprimary:
            l = ('primary', mass)
            m.write(l)
        if rewriteit:
            m.rewrite()
        print(mol + ': F_cal = ' + str(F_cal))
    
    if trendline:
        rsf_max = max([m.rsf for mol, m in mdict.items()])
        ax.plot([0, rsf_max], [0, r * rsf_max], 'k--')
        
    ax.set_xlabel('Relative Sensitivity Factor / [' + RSF_unit +']')
    ax.set_ylabel('F_cal / [C/mol]')
    #tickmol = ['C2H6', 'C2H4', 'CH4', 'H2', 'CO2']
    #ticknr = [0, 1, 1.5]
    #ax.set_xticks(ticknr + [RSF_dict[i] for i in tickmol])
    #ax.set_xticklabels([str(n) for n in ticknr] +tickmol)
    #ax.set_ylim([0,20])

    print('\nfunction \'recalibrate\' finished!\n\n')
    return mdict, ax



def get_signal(MS_data, mass, tspan='tspan_2', removebackground=False, 
             unit='A', verbose=True):
    '''
    Returns [x, y] where x is the time and y is QMS signal.
    A bit trivial, but I like having this function to work in parrallel 
    to get_flux.
    '''    

    if verbose:
        print('geting signal for ' + mass)
    
    x = MS_data[mass + '-x']
    y = MS_data[mass + '-y']

    if unit[-1] == 'A':
        if unit[:-1] == 'n' or unit[:-1] == 'nano':
            y = y*1e9
        elif unit[:-1] == 'u' or unit[:-1] == 'micro':
            y = y*1e6
    
    if tspan is None:
        tspan = 'tspan'
    if type(tspan) is str and not tspan=='all':
        tspan = MS_data[tspan]
    if not tspan == 'all':
        x, y = cut(x,y,tspan)  
    
    if removebackground:
        if type(removebackground) is float:
            background = removebackground * min(y)
        else:
            background = min(y)
        y = y - 0.99*background #Removing the entire background would fuck with log scales.
    return [x, y]


def get_flux(MS_data, mol, tspan='tspan_2',  
             unit='pmol/s', verbose=True, 
             removebackground=False, background='constant', endpoints=3):
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
    
    if tspan is None:
        tspan = 'tspan'
    if type(tspan) is str and not tspan=='all':
        tspan = MS_data[tspan]
    if not tspan == 'all':
        x, y = cut(x,y,tspan) 
    
    if removebackground:
        if background=='constant':
            if type(removebackground) is float:
                background = removebackground * min(y)
            else:
                background = min(y)
        elif type(background) is float:
            #background = background
            pass
        elif background=='linear':
            x_end = [np.average(x[:endpoints]), np.average(x[-endpoints:])]
            y_end = [np.average(y[:endpoints]), np.average(y[-endpoints:])]
            background = np.interp(x, x_end, y_end) 
        
        y = y - 0.99*background #so that we don't break the log scale.
        #I should get rid of this and assume the caller knows what they're doing.
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
    
    from EC_MS import set_figparams
    
    set_figparams(figwidth=8)
    
    F_cals = {'CO2': 23.371483271036237, 'H2': 29.899494846007542, 'O2': 14.756572997297784}
    mdict, ax = RSF_to_F_cal(calmol = {'CO2':'M44'}, RSF_source='Hiden', 
                             trust=['H2','O2','CO2'], F_cals=F_cals)
    #calmol is the molecule(s) who's calibration is extrapolated to the others
    #RSF_source can be 'NIST' or 'Hiden'.
    #trust gives the molecules that become the red dots.
    #F_cals overrides the stored calibration factors for given molecules.
    #mdict containes the Molecule objects with the extrapolated F_cal
    #ax is the axis it's plotted on.
    
    for mol, m in mdict.items():
        rsf = m.get_RSF(RSF_source='Hiden', verbose=False)
        print(mol + '\tF_cal = ' + str(m.F_cal) + '\tRSF = ' + str(rsf))
     
    ax.set_xlabel('Relative sensitivity factor / [a.u.]')
    ax.set_ylabel('F$_{\mathrm{cal}}$ / [C mol$^{-1}$]')
    
    plt.savefig('Fcal_vs_RSF.png')
    
        
    
    
    
    