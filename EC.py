

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
#import os

from .Data_Importing import epoch_time_to_timestamp
from .Combining import cut_dataset, is_EC_data, get_timecol, is_time, get_type


E_string_list = ['Ewe/V', '<Ewe>/V', '|Ewe|/V']
V_string_default = 'U vs RHE / [V]'
I_string_list = ['I/mA', '<I>/mA', '|EI|/mA']
J_string_default = 'J / [mA/cm^2]'


EC_cols_0 = ['mode', 'ox/red', 'error', 'control changes', 'time/s', 'control/V',
           'Ewe/V', '<I>/mA', '(Q-Qo)/C', 'P/W', 'loop number', 'I/mA', 'control/mA',
           'Ns changes', 'counter inc.', 'cycle number', 'Ns', '(Q-Qo)/mA.h',
           'dQ/C', 'Q charge/discharge/mA.h', 'half cycle', 'Capacitance charge/µF',
           'Capacitance discharge/µF', 'dq/mA.h', 'Q discharge/mA.h', 'Q charge/mA.h',
           'Capacity/mA.h', 'file number', 'file_number', 'Ece/V',
           'Ewe-Ece/V', '<Ece>/V', '<Ewe>/V', 'Energy charge/W.h', 'Energy discharge/W.h',
           'Efficiency/%', 'Rcmp/Ohm', 'time/s*', 'selector', 'j / [A/mg]',
           V_string_default, J_string_default,
           'I/ECSA / [mA/m^2]'
           ] # exotic J_str's. I need to change how this works!]


def select_cycles(EC_data_0, cycles=1, t_zero=None, verbose=True,
                  cycle_str=None, cutMS=True, data_type='CV', override=False):
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
    EC_data['data_cols'] = EC_data['data_cols'].copy()
    # ^ otherwise this gives me a problem elsewehre

    #it looks like I actually want Ns for CA's and cycle number for CV's.
    #How to determine which
    if cycle_str is None:
        if 'cycle number' in EC_data['data_cols']:
            cycle_str = 'cycle number'
        elif 'Ns' in EC_data['data_cols']:
            cycle_str = 'Ns'
        else:
            print('no cycle numbers detected!')

    cycle_numbers = EC_data[cycle_str]


    #N = len(cycle_numbers) #unneeded
    if type(cycles)==int:
        cycles = [cycles]
    mask = np.any(np.array([cycle_numbers == c for c in cycles]), axis=0)
    #list comprehension is awesome.

    for col in EC_data['data_cols']:
        try:
            if not get_type(col) in ['MS','cinfdata']:
                #then we're dealing with EC data
                try:
                    EC_data[col] = EC_data[col].copy()[mask]
                except KeyError:
                    print('hm... \'' + col + '\' in EC_data[data_cols] but not in EC_data')
        except IndexError:
            print('trouble selecting cycle ' + str(cycles) + ' of ' + col + '\n' +
                    'type(mask) = ' + str(type(mask)))
            good = False
    t0 = EC_data['tstamp']
    tspan = np.array([min(EC_data['time/s']), max(EC_data['time/s'])])

    if cutMS:
        time_masks = {}
        for col in EC_data['data_cols']:
            if not is_EC_data(col): #then we've got a QMS or synchrotron variable
                timecol = get_timecol(col)
                print(col + '  ' + timecol) #debugging
                if timecol in time_masks:
                    mask = time_masks[timecol]
                else:
                    try:
                        mask = np.logical_and(tspan[0] < EC_data[timecol], EC_data[timecol] < tspan[-1])
                        time_masks[timecol] = mask
                    except KeyError:
                        print('can\'t cut ' + col + ', don\'t have timecol ' + timecol + '. skipping.')
                        continue
                #if not is_time(col):
                EC_data[col] = EC_data[col][mask]

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

        tspan = tspan - t_zero
        t0 = t0 + t_zero
        EC_data['tstamp'] = t0
    EC_data['timestamp'] = epoch_time_to_timestamp(t0)
    EC_data['tspan'] = tspan
    EC_data['tspan_2'] = tspan
    EC_data['tspan_0'] = tspan + t0
    EC_data['data_type'] += ' selected'
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
        # Colors in the area between two cycles in the specified
        # potential range (Vspan) and direction (redox=1 for
        # anodic.)
        # Returns (dQ, data), where dQ is the difference in
        # charge passed during that region, and data=[t, V, J_diff]
        # has columns for the time, potential, and difference in
        # current for the specified region.
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
        #print(type(cycles_data)) # debugging
        V_str, J_str = sync_metadata(cycle_data, verbose=verbose)
        print(V_str + ', ' + J_str) # debugging
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
            print('didn''t find A_el. Using A_el=1')
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
                t_i=0, redox_str='ox/red', verbose=True, closecycle=False):
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

    I_start = np.argmax(t>t_i) #so that I get point 0 in the first cycle.
    I_finish = I_start + 1
    I_next = I_start + 1
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


    try:
        return [cyclesets[i] for i in cycles] #Whoa.
    except:
        return cyclesets

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



def make_selector(data, sel_str='selector'):
    changes = np.tile(False, data['time/s'].shape)
    for col in ['cycle number', 'loop number', 'file number']:
        if col in data:
            n = data[col]
            n_up = np.append(n[1:], n[-1])
            changes = np.logical_or(changes, n_up>n)
    selector = np.cumsum(changes)
    data[sel_str] = selector
    return sel_str


def sync_metadata(data, RE_vs_RHE=None, A_el=None,
                  V_str=None, J_str=None, E_str=None, I_str=None,
                  verbose=True):
    '''
    A nice one-serve-all function for calibrating data, updating calibration, checking
    if it's calibrated, and checking what the most useful data available is called.
    It does simple things, but it does a lot.
    Intuitiveness of use is valued over intuitiveness of internal working.
    Here's the details:

    -----

    data is a dictionary containing uncalibrated and perhaps calibrated data, as well
    as the name of the keys pointing to certain types of data.

    E_str and I_str are the keys pointing to uncalibrated variable 1 (i.e.
    potential vs RE) and variable 2 (i.e. absolute current), respectively.

    V_str and J_str are the keys pointing to the calibrated variable 1 (i.e.
    potential vs RHE) and variable 2 (i.e. normalized current), respectively.

    RE_vs_RHE, if given, calibrates variable 1 such that data[V_str] = data[E_str] + RE_vs_RHE

    A_el, if given, calibrates variable 2 such that data[J_str] = data[I_str] / A_el

    RE_vs_RHE and A_el are also stored in data as data['RE_vs_RHE'] and data['A_el'].
    If either RE_vs_RHE or A_el are set here to 'existing', the data is re-calibrated
    according to the stored value.

    All strings E_str, V_str, I_str, and J_str are pointed to by data['E_str'] etc after
    a call to this function.

    All strings E_str, V_str, I_str, and J_str can be updated to point to a new column in the dataset.
    This does not recalibrate the data unless the relavent of RE_vs_RHE or A_el are given.

    If RE_vs_RHE and/or A_el are given and V_str and J_str are neither specified here
    nor pre-specified in data['V_str'] or data['J_str'], new columns are created
    with the default names V_str = 'U vs RHE / [V]' and J_str = 'J / [mA/cm^2]'
    to contain the calibrated values for variable 1 and varible 2, respectively

    V_str are J_str are returned if known.
    If V_str is not known (no variable 1 calibration), E_str is returned in its place.
    If J_str is not known (no variable 2 calibration), I_str is returned in its place.

    if verbose=True, the function talks to you.

    '''

    if verbose:
        print('\nsyncing metadata for ' + data['title'] + '\n')

    # use these to keep track of whether variables 1 and 2 will be calibrated:
    cal1 = False
    cal2 = False

    # obviously we will calibrate if given RE_vs_RHE and/or A_el
    if RE_vs_RHE is not None:
        data['RE_vs_RHE'] = RE_vs_RHE
        cal1 = True
    if A_el is not None:
        data['A_el'] = A_el
        cal2 = True

    # Now, look through the strins input here and in data. Are we setting or
    # getting column names?

    if E_str is None and 'E_str' in data:   # then get it
        E_str = data['E_str']
    elif E_str in data:                     # then set it
        data['E_str'] = E_str

    if I_str is None and 'I_str' in data:   # then get it
        I_str = data['I_str']
    elif I_str in data:                     # then set it
        data['I_str'] = I_str

    if V_str is None and 'V_str' in data:   # then get it
        V_str = data['V_str']
    elif V_str in data:                     # then set it
        data['V_str'] = V_str
    elif V_str is not None:                 # A brand new V_str demands calibration of variable 1
        cal1 = True

    if J_str is None and 'J_str' in data:   # then get it
        J_str = data['J_str']
    elif J_str in data:                     # then set it
        data['J_str'] = J_str
    elif J_str is not None:                 # A brand new J_str demands calibration of variable 1
        cal2 = True

    # If we're calibrating variable 1 or still looking for a V_str, we need an E_str
    if E_str is None and (cal1 or V_str is None):   # then we'll need an E_str!
        try:   # see if a possible E_str (listed at the top of the module) is represented in the data
            E_str = next(s for s in E_string_list if s in data['data_cols'])
        except StopIteration:
            print('sync metadata can\'t find any value for E_str!')
            print('if you needed to calibrate variable 1, that won\'t happen now.')
            cal1 = False

    # If we're calibrating variable 2 or still looking for a J_str, we need an I_str
    if I_str is None and (cal2 or J_str is None):   # then we'll need an E_str!
        try:   # see if a possible I_str (listed at the top of the module) is represented in the data
            I_str = next(s for s in I_string_list if s in data['data_cols'])
        except StopIteration:
            print('sync metadata can\'t find any value for I_str!')
            print('if you needed to calibrate variable 2, that won\'t happen now.')
            cal2 = False

    #----------  alright, now we're ready to calibrate! --------------

    if cal1:   #Calibrate variable 1 (i.e., electrode potential)
        E = data[E_str]            # get the data to be calibrated
        # get the calibration factor
        if (RE_vs_RHE is None or RE_vs_RHE =='existing') and 'RE_vs_RHE ' in data:
            A_el = data['A_el']
        if RE_vs_RHE  is None:
            print('Your call to sync_metadata demands calibration of variable 2, ' +
                  'but I can\'t figure out what RE_vs_RHE is. Using RE_vs_RHE = 0.')
            RE_vs_RHE = 0
        V = E + RE_vs_RHE          #calibrate the data
        if V_str is None:         # figure out where to put it
            V_str = V_string_default
        data['V_str'] = V_str   # remind yourself where you're putting it
        data[V_str] = V         # and put it there!


    if cal2:  #Calibrate variable 2 (i.e., electrical current)
        I = data[I_str]              # get the data to be calibrated
        # get the calibration factor
        if (A_el is None or A_el=='existing') and 'A_el' in data:
            A_el = data['A_el']
        if A_el is None:
            print('Your call to sync_metadata demands calibration of variable 2, ' +
                  'but I can\'t figure out what A_el is. Using A_el=1.')
            A_el = 1
        J = I / A_el        #calibrate the data
        if J_str is None:        # figure out where to put it
            J_str = J_string_default
        data['J_str'] = J_str  # remind yourself where you're putting it
        data[J_str] = J        # and put it there!

    #----------  and, let's get out of here! --------------

    # make sure this function points to uncalibrated data if necessary.
    data['E_str'] = E_str
    data['I_str'] = I_str
    if V_str is None:
        V_str = E_str
    if J_str is None:
        J_str = I_str

    # oh, yeah, and they're all data columns!
    for s in [E_str, I_str, V_str, J_str]:
        if s is not None and s not in data['data_cols']:
            data['data_cols'] += [s]


    # return the keys for the most useful data
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





def get_through_sweep(data=None, t_str=None, V_str=None, t=None, V=None, t_i=0,
                      Vspan=[0.4, 0.6], edge=0.01, redox=None, verbose=True):
    '''
    Return i_start and i_finish, the indeces corresponding to the first
    complete anodic(redox=True) or cathodic (redox=False) sweep through V_span
    starting after t_i. t and V can be given directly, or
    V_str and J_str can point to the corresponding columns in data.
    '''
    # parse inputs:
    if t is None:
        if t_str is None:
            t_str = data['t_str']
        t = data[t_str]
    if V is None:
        if V_str is None:
            V_str = data['V_str']
        V = data[V_str]
    if redox is None: #then assume they gave it in the order of Vspan.
        redox = Vspan[0] < Vspan[-1]
    #print(redox) #debugging

    # define some things to generalize between anodic and cathodic
    def before(a, b):
        if redox:
        # before means more cathodic if we want the anodic sweep
            return a < b
        else:
        # and more anodic if we want the cathodic sweep
            return a > b
    if redox:
        # we start with the lower limit of V_span if we want the anodic sweep
        Vspan = np.sort(np.array(Vspan))
        Vspan_wide = [Vspan[0]-edge, Vspan[-1]+edge]
    else:
        # and with the upper limit of V_span if we want the cathodic sweep
        Vspan = - np.sort(-np.array(Vspan))
        Vspan_wide = [Vspan[0]+edge, Vspan[-1]-edge]

    #print('len(t) = ' + str(len(t))) # debugging
    #print('len(V) = ' + str(len(V))) # debugging
    #print(str(len(t>t_i)) + ' ' + str(len(before(Vspan[0], V)))) #debugging

    t_out = t[np.argmax(np.logical_and(t > t_i,
                                       before(V, Vspan_wide[0])
                                       ) # True if after t_i and comfortably out on start side
                        ) #first index for which V is comfortably out on start side
             ] #corresponding time
    #print(f't_i = {t_i}, t_out = {t_out}, t={t}') # debugging
    i_start = np.argmax(np.logical_and(t>t_out, before(Vspan[0], V)))
    # ^ first index of full sweep through range
    i_finish = np.argmax(np.logical_and(t>t_out, before(Vspan[1], V))) - 1
    # ^ last index of full sweep through range
    if verbose:
        print('first scan from ' + str(Vspan[0]) + ' --> ' +
              str(Vspan[1]) + ' occurs between indeces ' + str(i_start) +
              ' and ' + str(i_finish))
    return i_start, i_finish


def get_shunt_current_line(data, V_DL, t_i=0,
                           t_str=None, V_str=None, I_str=None, N=100, ax=None,
                           verbose=True):
    '''

    '''
    if t_str is None:
        try:
            t_str = data['t_str']
        except KeyError:
            t_str = 'time/s'
    if V_str is None:
        V_str = data['V_str']
    if I_str is None:
        I_str = data['I_str'] # I don't want to deal with area normalization.

    if 'RE_vs_RHE' in data:
        V_unit = 'V vs RHE'
    else:
        V_unit = 'V'

    V = data[V_str]
    I = data[I_str]
    t = data[t_str]
    #print('t_str = ' + t_str + ', V_str = ' + V_str) # debugging

    an_start, an_finish = get_through_sweep(t=t, V=V, Vspan=V_DL, redox=True, t_i=t_i)
    cat_start, cat_finish = get_through_sweep(t=t, V=V, Vspan=V_DL, redox=False, t_i=t_i)

    V_an = V[an_start:an_finish]
    I_an = I[an_start:an_finish]

    V_cat = V[cat_start:cat_finish]
    I_cat = I[cat_start:cat_finish]

    Vspan = [max(min(V_an), min(V_cat)), min(max(V_an), max(V_cat))]
    V_interp = np.arange(N)/N * (Vspan[-1] - Vspan[0]) + Vspan[0]

    I_an_interp = np.interp(V_interp, V_an, I_an)
    fixcat = np.argsort(V_cat)
    I_cat_interp = np.interp(V_interp, V_cat[fixcat], I_cat[fixcat])

    I_avg = (I_an_interp + I_cat_interp) / 2

    pfit = np.polyfit(I_avg, V_interp, deg=1)

    t_f = t[max(an_finish, cat_finish) + 1]

    if verbose:
        print('function \'get_shunt_current_line\' calculates R_shunt = ' +
              str(pfit[0]) + ' kOhm, with potential of zero shunt = ' +
              str(pfit[1]) + ' ' + V_unit)
    if ax is not None:
        try:
            A_el = data['A_el']
            J_str = data['J_str']
        except KeyError:
            A_el = 1
            J_str = I_str
        if ax == 'new':
            fig, ax = plt.subplots()
        ax.plot(V_interp, I_avg/A_el, 'r')
        ax.plot(V_an, I_an/A_el, 'k')
        ax.plot(V_cat, I_cat/A_el, 'k')
        ax.plot(V_interp, I_an_interp/A_el, 'k--')
        ax.plot(V_interp, I_cat_interp/A_el, 'k--')
        x = [np.min(V_interp), np.max(V_interp)]
        y = (x - pfit[1])/pfit[0]
        ax.plot(x, y/A_el, 'g')
        ax.set_xlabel(V_str)
        ax.set_ylabel(J_str)
    return pfit, t_f



def correct_shunt(data, tspan='all', R_shunt=None, V_intercept=None, pfit=None,
                  t_str=None, V_str=None, I_str=None, J_str=None, verbose=True,
                  **kwargs):
    '''
    '''
    if verbose:
        print('correcting shunt!')
    if pfit is None and R_shunt is None and V_intercept is None:
        #print(I_str)  #debugging
        pfit, t_f = get_shunt_current_line(data,
                                      t_str=t_str, V_str=V_str, I_str=I_str,
                                      verbose=verbose, **kwargs)
    if R_shunt is None:
        R_shunt = pfit[0]  # CE-ground shunt resistance in kOhm
    if V_intercept is None:
        V_intercept = pfit[1] # potential vs RHE of zero shunt current

    if t_str is None:
        try:
            t_str = data['t_str']
        except KeyError:
            t_str = 'time/s'
    if V_str is None:
        V_str = data['V_str']
    if I_str is None:
        I_str = data['I_str']
    if J_str is None:
        J_str = data['J_str']

    t = data[t_str]
    V = data[V_str]
    I = data[I_str].copy() #preserve the un-corrected originals
    J = data[J_str].copy()
    A_el = data['A_el']

    t = data[t_str]
    if tspan is None or tspan=='all':
        mask = np.tile(True, np.shape(t))
    else:
        mask = np.logical_and(tspan[0]<t, t<tspan[-1])

    I_shunt = 1/R_shunt * (V[mask] - V_intercept)
    I[mask] = I[mask] - I_shunt
    J[mask] = J[mask] - I_shunt/A_el

    I_str += '*'
    J_str += '*'
    if verbose:
        print('saving corrected data as ' + J_str + ' (normalized current) and ' +
              I_str + ' (raw current) + ')
    data[I_str] = I
    data[J_str] = J
    data['I_str'] = I_str
    data['J_str'] = J_str
    if I_str not in data['data_cols']:
        data['data_cols'] += [I_str]
    if J_str not in data['data_cols']:
        data['data_cols'] += [J_str]

    return pfit


def get_scan_rate(data, Vspan=[0.3, 0.6], V_str=None, J_str=None, t_i=0):
    '''
    Returns scan rate in mV/s. Needs a Vspan.
    '''
    if V_str is None or J_str is None:
        V_str_1, J_str_1 = sync_metadata(data)
    if V_str is None:
        V_str = V_str_1
    t_str = 'time/s'

    t, V = data[t_str], data[V_str]

    I_start_an, I_finish_an = get_through_sweep(data, Vspan=Vspan, t_i=t_i, redox=1)

    scan_rate = (V[I_finish_an] - V[I_start_an])/(t[I_finish_an] - t[I_start_an])*1e3 # mV/s

    return scan_rate

def get_capacitance(data, V_DL=[0.3, 0.6], V_str=None, J_str=None, t_i=0):
    '''
    Returns capacitance in F/cm^2
    '''
    if V_str is None or J_str is None:
        V_str_1, J_str_1 = sync_metadata(data)
    if V_str is None:
        V_str = V_str_1
    if J_str is None:
        J_str = J_str_1
    t_str = 'time/s'

    t, V, J = data[t_str], data[V_str], data[J_str]

    I_start_an, I_finish_an = get_through_sweep(data, Vspan=V_DL, t_i=t_i, redox=1)
    J_an = np.mean(J[I_start_an:I_finish_an]) # average current on anodic sweep in mA/cm^2

    I_start_cat, I_finish_cat = get_through_sweep(data, Vspan=V_DL, t_i=t_i, redox=0)
    J_cat = np.mean(J[I_start_cat:I_finish_cat]) # average current on cathodic sweep in mA/cm^2

    scan_rate = (V[I_finish_an] - V[I_start_an])/(t[I_finish_an] - t[I_start_an])*1e3 # mV/s

    cap = (J_an - J_cat)/2 / scan_rate # [mA/cm^2] / [mV/s] = [ (C/s)/(V/s) /cm^2] = F/cm^2
    return cap

def correct_ohmic_drop(data, R_ohm=0):
    V_str, J_str = sync_metadata(data)
    I_str = data['I_str']
    V, I = data[V_str].copy(), data[I_str].copy()
    V -= R_ohm * I*1e-3
    if not V_str[-1] == '*':
        V_str = V_str + '*'
    data[V_str] = V
    data['V_str'] = V_str
    return V_str