

# -*- coding: utf-8 -*-
"""
created 16I15
last edited 16J27

@author: Scott

functions for use on EC-lab data
"""

#make python2-compatible
from __future__ import print_function
from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
import os

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from Combining import timestamp_to_seconds, is_time, cut 
    from Combining import cut_dataset#, get_time_col
else:                           #then we use relative import
    from .Combining import timestamp_to_seconds, is_time, cut 
    from .Combining import cut_dataset#, get_time_col


def select_cycles(EC_data_0, cycles=1, t_zero=None, verbose=True, cycle_str=None, cutMS=True, data_type='CV', override=False):
    ''' 
    This function selects one or more cycles from EC_data_0.
    Use this before synchronizing!
    Works for both CA and CV
    #changed 16L22 to work on EC_and_MS data
    #just set cycle_str to 'loop number' to select loop rather than cycle.
    #override to ignore when cut returns empty dataset.
    '''
    if verbose:
        print('\nSelecting cycles ' + str(cycles) + ' from \'' + EC_data_0['title'] + '\'\n')
    
    good = True    
    EC_data = EC_data_0.copy()
    
    #it looks like I actually want Ns for CA's and cycle number for CV's.
    #How to determine which
    if cycle_str is None:
        if data_type == 'CV' and 'cycle number' in EC_data['data_cols']:
            cycle_str = 'cycle number'
        elif 'Ns' in EC_data['data_cols']:
            cycle_str = 'Ns'
        else:
            print('no cycle numbers detected!')

    cycle_numbers = EC_data[cycle_str]


    N = len(cycle_numbers)
    if type(cycles)==int:
        cycles = [cycles]
    I_keep = np.array([I for I in range(N) if cycle_numbers[I] in cycles])
    #list comprehension is awesome.

    for col in EC_data['data_cols']:
        try:
            if not (col[0] == 'M' and col[-2:] in ['-x', '-y']):  
                #then we're dealing with EC data
                try:
                    EC_data[col] = EC_data[col].copy()[I_keep]
                except KeyError:
                    print('hm... \'' + col + '\' in EC_data[data_cols] but not in EC_data')
        except IndexError:
            print('trouble selecting cycle ' + str(cycles) + ' of ' + col + '\n' +
                    'type(I_keep) = ' + str(type(I_keep)))
            good = False
    t0 = timestamp_to_seconds(EC_data['timestamp'])
    tspan = np.array([min(EC_data['time/s']), max(EC_data['time/s'])])
    EC_data['tspan'] = tspan
    EC_data['tspan_2'] = tspan
    EC_data['tspan_0'] = tspan + t0
    EC_data['data_type'] += ' selected'   
    
    if cutMS:
        for col in EC_data['data_cols']:
            #print(col)
            if col[0] == 'M' and col[-2:] == '-x': #then we've got a QMS time variable
                y_col = col[:-2] + '-y'
                #print('select cycles is cutting in MS data ' + col )
                EC_data[col], EC_data[y_col] = cut(EC_data[col], EC_data[y_col], tspan, override=override)       

    if t_zero is not None:
        if verbose:
            print('\'select_cycles\' is resetting t_zero')
            print('t_zero = ' + str(t_zero))
        if type(t_zero) is str:
            try:
                n = eval(t_zero)
                #e.g. t_zero = '3' sets t=0 to the start of the third cycle, 
                #regardless of the selected cycles 
                if type(n) is not int:
                    raise NameError
                t_zero = next(EC_data['time/s'][i] for i,c in enumerate(EC_data[cycle_str]) if c==n)

            except NameError:     
                #this should be the case if t_zero=='start'
                t_zero = tspan[0]
            if verbose:
                print('aka, shifting by t_zero=' + str(t_zero))
                
        for col in EC_data['data_cols']:
            if is_time(col):
                EC_data[col] = EC_data[col] - t_zero

        EC_data['tspan'] = tspan - t_zero #fixed from tspan - tspan[0] - t_zero 17H09
        EC_data['tspan_2'] = EC_data['tspan']
        
    EC_data['good'] = good
    return EC_data


def remove_delay(CV_data):
    '''
    Gets rid of the delay at the beginning of .mpt files before it actually starts 
    cycling. This is not seen cycle_number, but in control changes, which goes to 0 for 
    the first time right as the cycle starts... I think. 16I29
    '''
    control = CV_data['control changes']
    
    I_start = np.where(control == 0)[0][0]

    for col in CV_data['data_cols']:
        CV_data[col] = CV_data[col][I_start:-1]

    return CV_data
        

def CV_difference(cycles_data, redox=1, Vspan=[0.5, 1.0], 
                  ax=None, color='g', verbose=True):
    '''
    This will calculate the difference in area between two cycles in a CV, 
    written for CO stripping 16J26. If ax is given, the difference will be
    filled in with color.
    '''
    if verbose:
        print('\n\nfunction \'CV_difference\' at your service!\n')  
    
    if redox == 'ox':
        redox = [1]
    elif redox == 'red':
        redox = [0]
    elif type(redox) is int:
        redox = [redox]
    elif redox is None:
        redox = [0,1]
    
    Vs = []   
    Js = []
    Q = []
    JV = []
    ts = []
    for cycle_data in cycles_data:
        #print(type(cycles_data))
        V_str, J_str = sync_metadata(cycle_data, verbose=verbose)
        V = cycle_data[V_str]
        J = cycle_data[J_str]
        t = cycle_data['time/s']
        
        ro = cycle_data['ox/red']
        q = cycle_data['(Q-Qo)/C']
        I_keep = [I for (I, (V, ro)) in enumerate(zip(V, ro)) if 
                    Vspan[0] < V < Vspan[1] and ro in redox]
        
        V = V[I_keep]        
        J = J[I_keep]
        t = t[I_keep]
        print('V_range starts at t = ' + str(t[0]))
        
        Vs += [V]
        Js += [J]
        ts += [t]
        
        Q += [q[I_keep[-1]] - q[I_keep[0]]]       
        JV += [np.trapz(J, V)]
        
    dQ = Q[0] - Q[1] 
    dJV = JV[0] - JV[1]
    
    if verbose:
        try:
            A_el = cycle_data['A_el']
        except KeyError:
            A_el = 1
            print('didn''t find A_el.')
        if A_el is None:
            A_el = 1
        print('difference in charge passed: a = ' + str(dQ) + ' C\n' + 
                'difference in CV area: b = ' + str(dJV) + ' V*mA/cm^2\n' + 
                'This implies a scan rate of: b/a*A_el = ' + str(dJV / dQ * A_el) + ' mV/s') 
    
    # We're going to return three vectors, for t V, and J, and 
    #   all of them will be the same length as the first dataset, i.e. t[0]
    if len(Vs[0]) != len(Vs[1]):  #then we'll have to interpolate
        if 1 in redox:
            Js[1] = np.interp(Vs[0], Vs[1], Js[1])
            V = Vs[0]
        else:    
            Js[1] = np.interp(-Vs[0], -Vs[1], Js[1])
            V = Vs[0]
    else:
        V = (Vs[0] + Vs[1]) / 2
    J_diff = Js[0] - Js[1] #note this is all optimized for CO stripping
    t = ts[0]
    
    if ax:
        if ax == 'new':
            ax = plt.figure().add_subplot(111)
            ax.set_xlabel(V_str)
            ax.set_ylabel(J_str)
        ax.fill_between(V, Js[0], Js[1], where=Js[0]>Js[1],
                        facecolor=color, interpolate=True)
    if verbose:
        print('\nfunction \'CV_difference\' finished!\n\n')
    
    return dQ, [t, V, J_diff]


def clip_cycles(dataset, cycles=1, V_clip=0, redox=1, V_str=None, t_str='time/s',
                redox_str='ox/red', verbose=True, closecycle=False):
    '''
    puts the clip at a specified potential (or other data column given 
    by V_str) V_clip, and returns a subset given by indeces in cycles. By
    default returns the first full cycle in the dataset.
    if redox=1, cuts on the anodic sweep, if redox=0 on the cathodic sweep.
    '''
    print(redox)
    if verbose:
        print('\n\nfunction \'clip_cycles\' at your service!\n')       
    
    if V_str is None:
        V_str, J_str = sync_metadata(dataset, verbose=False)
        
    if type(cycles) is int:  #my need to always do this kind of thing is 
        cycles = [cycles]    #an annoying aspect of python.
    
    t, V, ro = dataset[t_str].copy(), dataset[V_str].copy(), dataset[redox_str].copy() 
    #wouldn't want these to get fucked up
    N = len(V)

    if redox: #I think this is more efficient than putting the if inside the
    #function, because it doesn't have to keep reevaluating truth value of redox
        print('I_finish will be when redox==1.')
        def condition(I):
            return V[I] > V_clip and ro[I] == 1
        #V[I+1] > V[I] doesn't always work. 
    else:
        print('I_finish will be when redox==0.')
        def condition(I):
            return V[I] < V_clip and ro[I] == 0
    n = 0

    I_start = 0 #so that I get point 0 in the first cycle.
    I_finish = 1
    I_next = 1
    cyclesets = []
    endit = False
    while n < max(cycles) + 1:
        print('I_start = ' + str(I_start))        
        print('t[I_start] = ' + str(t[I_start]))
        print('V[I_start] = ' + str(V[I_start]))
        
        I_next = next(I for I in range(I_start+1, N) if not condition(I))
        print('I_next = ' + str(I_next))        
        print('t[I_next] = ' + str(t[I_next]))
        print('V[I_next] = ' + str(V[I_next]))
        #Choose I_next to be on the subsequent scan, so that I don't just
        #cut it into a lot of single points.
        
        try:
            I_finish = next(I for I in range(I_next, N) if condition(I))
            print('I_finish = ' + str(I_finish))        
            print('t[I_finish] = ' + str(t[I_finish]))
            print('V[I_finish] = ' + str(V[I_finish]))
        except StopIteration:
            print('StopIteration')
            I_finish = N-1
            endit = True
        except IndexError:
            print('IndexError')
            endit = True
            I_finish = N-1        
        tspan = [t[I_start], t[I_finish]]
        if not tspan[1] > tspan[0]:
            print('warning! tspan = ' + str(tspan))
        print('cutting dataset')
        c = cut_dataset(dataset, tspan)
        if closecycle:
            c = closecycle[c]
        cyclesets += [c]
        print('got a cycle! len(cyclesets) = ' + str(len(cyclesets)) + '\n')
        if endit:
            print('but also hit a problem. We\'re done here.')
            break
        I_start = I_finish
        n += 1
        print('\n\n')

        
    if len(cycles) == 1:
        try:
            return cyclesets[cycles[0]]
        except IndexError:
            print('couldn\'t get your cycle. returning the first one.')
            return cyclesets[0]

    if verbose:
        print('\nfunction \'clip_cycles\' finished!\n\n')       
 


    return [cyclesets[i] for i in cycles] #Whoa.

def close_cycle(cycle_0):
    '''
    joins the ends of the data in a cycle to make it 
    look nice when plotted vs potential.
    '''
    cycle = cycle_0.copy()
    for col in cycle['data_cols']:
        x = cycle[col]
        if is_time(col):
            x = np.append(x,2 * x[-1] - x[-2]) #continue same t spacing
        else:
            x = np.append(x, x[0])
        cycle[col] = x
    return cycle
    
def smooth_pulses(CA_Data_0, verbose=1):
    '''
    This function turns the CA data into a square wave by averaging the
    potential over the duration of a pulse (where it should be constant).
    Useful when noise makes the figures look ugly dispite otherwise good data.
    If you have to use this function, though, I would say the results are not 
    publication-ready.
    '''
    if verbose:
        print('\n\nfunction \'smooth_pulses\' at your service!')
    CA_Data = CA_Data_0.copy()
    cycle_numbers = CA_Data['Ns']
    cycles = np.unique(cycle_numbers)
    for c in cycles:
        I_cycle = np.array([i for (i,cycle) in enumerate(cycle_numbers) if cycle==c])
        V_avg = np.average(CA_Data['Ewe/V'][I_cycle])
        CA_Data['Ewe/V'][I_cycle] = V_avg
    if verbose:
        print('function \'smooth_pulses\' finished!\n\n')
    return CA_Data


def plot_vs_time(EC_data, axes='new', y_strings='default', colors=None,
                 RE_vs_RHE=None, A_el=None, verbose=1):
    
    if verbose:
        print('\n\nfunction \'plot_vs_time\' at your service!')    
    
    V_str, J_str = sync_metadata(EC_data, RE_vs_RHE, A_el)
    if y_strings == 'default':
        y_strings = V_str, J_str
    
    t_str = 'time/s'
    t = EC_data[t_str]
    
    if colors is None:
        colors = ['k'] * len(y_strings)
    
    if axes == 'new':
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax2 = ax1.twinx()
        axes = [ax1, ax2]
    
    for (ax, y_str, color) in zip(axes, y_strings, colors):
        try:
            y = EC_data[y_str]
        except KeyError:
            print('Can''t find ' + y_str + '. skipping that one.')
            continue
        
        ax.plot(t, y, color, label=y_str)
        ax.set_xlabel(t_str)        
        ax.set_ylabel(y_str)
    
    if verbose:
        print('function \'plot_vs_time\' finished!\n\n')  
        
    return axes


def sync_metadata(EC_data, RE_vs_RHE=None, A_el=None, verbose=True):
    '''
    Deal with all the annoying RE and J vs I vs <I> stuff once and for all here.
    After this has been called, all plotting methods need only to utilize
    EC_data['V_str'] and EC_data['J_str']
    '''    
    if verbose:
        print('\nsyncing metadata for ' + EC_data['title'] + '\n')
    
    if RE_vs_RHE is None and A_el is None and 'J_str' in EC_data.keys() and 'V_str' in EC_data.keys():
        #added 17G26 so that plot_experiment would stop rewriting data[V_str] and data[J_str]
        if verbose:
            print('... already sync\'d! \n\n')
        return EC_data['V_str'], EC_data['J_str'] 
    
    if RE_vs_RHE is not None:
        EC_data['RE_vs_RHE'] = RE_vs_RHE
    elif 'RE_vs_RHE' in EC_data:
        RE_vs_RHE = EC_data['RE_vs_RHE']
    else:
        EC_data['RE_vs_RHE'] = None
    try:
        E_str = [s for s in ['Ewe/V', '<Ewe>/V', '|Ewe|/V'] if s in EC_data['data_cols']][0] 
    except IndexError:
        print('data doesn\'t include Ewe!')  
        E_str = None
        V_str = None
    if RE_vs_RHE is None:
        V_str = E_str
    elif E_str is not None:
        V_str = 'U vs RHE / [V]' #changed from E to U 17E21
        EC_data[V_str] = EC_data[E_str] + RE_vs_RHE
    
    if A_el is not None:
        EC_data['A_el'] = A_el
    elif 'A_el' in EC_data:
        A_el = EC_data['A_el']
    else:
        EC_data['A_el'] = None

    try:
        I_str = [s for s in ['I/mA', '<I>/mA', '|EI|/mA'] if s in EC_data['data_cols']][0] 
    except IndexError:
        print('data doesn\'t include I!')
        I_str = None
        J_str = None
    if A_el is None:
        J_str = I_str
    elif I_str is not None:
        J_str = 'J /[mA/cm^2]'
        EC_data[J_str] = EC_data[I_str] / A_el

    EC_data['E_str'] = E_str 
    EC_data['V_str'] = V_str   
    EC_data['J_str'] = J_str
    EC_data['I_str'] = I_str
    
    EC_data['data_cols'] = EC_data['data_cols'].copy() #17B02
    for col in [V_str, J_str]:
        if col not in EC_data['data_cols'] and col is not None:
            EC_data['data_cols'] += [col]
            if verbose:
                print('added ' + col + ' to data_cols')
        
    return V_str, J_str
    

def plot_CV_cycles(CV_data, cycles=[0], RE_vs_RHE=None, A_el=None, ax='new',
                   cycle_str='cycle number',
                   saveit=0, title='default', leg=0, verbose=1, colors=None):
    '''
    plots a subset of cycles in a CV
    '''
    if verbose:
        print('\n\nfunction \'plot_CV_cycles\' at your service!')
    
    if ax == 'new':
        fig1 = plt.figure()
        ax = fig1.add_subplot(111) 
    
    if type(cycles)==int:
        cycles = [cycles]  
        
    V_str, J_str = sync_metadata(CV_data, RE_vs_RHE, A_el) #added 16J26
    
    data_to_return = [] 
    for n, cycle in enumerate(cycles):

        cycle_data = select_cycles(CV_data, cycles=cycle, verbose=verbose, cycle_str=cycle_str)
        data_to_return += [cycle_data]  #added 16J25

        if ax is not None:
            V = cycle_data[V_str]
            J = cycle_data[J_str]
    
            if colors:
                color = colors[n]
            else: color = 'k'
            
            ax.plot(V, J, color) 
            
    if ax is not None:
        ax.set_xlabel(V_str)
        ax.set_ylabel(J_str)      
    
        if saveit:
            if title == 'default':
                title == CV_data['title'] + '.png'
            fig1.savefig(title)
            
    if verbose:
        print('\nfunction \'plot_CV_cycles\' finished!\n\n')  
    
    return data_to_return, ax


if __name__ == '__main__':
    from Data_Importing import import_data

    plt.close('all')    
    
    import_raw_data = 0
    if import_raw_data:
        EC_directory = ('/home/soren/Dropbox (Soren Scott Inc)' +
            '/Sniffer_Experiments/03_Pt_Sputtered/Data/16I28_for_Hawaii/')
        EC_file =  '03_HER_OER_C01.mpt'
        
        EC_data_0 = import_data(EC_directory + EC_file, data_type='EC')
    
    sync_metadata(EC_data_0, RE_vs_RHE=0.535, A_el=0.196)
    
    EC_data_1 = select_cycles(EC_data_0, [1,2], tzero='start')
    [ax1, ax2] = plot_vs_time(EC_data_1, colors=['r','k'])
    
    EC_data_2 = select_cycles(EC_data_0, [3,4], tzero='start')
    plot_vs_time(EC_data_2, colors=['r--','k--'], axes=[ax1, ax2])
    
    
    