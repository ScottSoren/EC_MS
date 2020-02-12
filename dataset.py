#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 18:12:03 2020

@author: scott
"""
import os, re, pickle
import numpy as np
from types import FunctionType
from functools import wraps
from matplotlib import pyplot as plt
from matplotlib import cm as colormap

from .EC import sync_metadata, make_selector, select_cycles
from .Data_Importing import load_from_file
from .Combining import synchronize, cut_dataset
from .Plotting import plot_experiment, plot_vs_potential
from .EC import correct_ohmic_drop, CV_difference


def get_data_from_file(file_name, data_type=None, verbose=True): # assumes you're already in the folder
    if type(file_name) is dict:
        return file_name # so that the dataset can be initiated with data already in a dictionary
    if re.search('.pkl$', file_name):
        with open(file_name, 'rb') as f:
            return pickle.load(f)
    elif re.search('.mpt$', file_name):
        return load_from_file(file_name, data_type='EC',
                                   verbose=verbose)
    elif data_type is not None:
        return load_from_file(file_name, data_type=data_type,
                                   verbose=verbose)
    else:
        print('WARNING: loading files of the type ' + file_name +
              ' is not yet implemented in Dataset.__init__() !!!'
              + ' Try specifying a data_type.'
              )

def with_update(method):
    @wraps(method)
    def method_with_update(self, *args, **kwargs):
        ret = method(self, *args, **kwargs)
        self.update_with_data()
        return ret
    return method_with_update


class Dataset:
    '''
    This class implements the dataset. Its design is to be back-compatible
    with the dataset dictionaries that were the main object in function-
    centric EC_MS programming.

    Dataset just serves as a wrapper around dataset dictionaries to
    make the package seem object-oriented. It has __getitem__ and
    __getattr__ methods that make (key,value) pairs and attributes somewhat
    interchangable. It should be back-compatable, but I haven't really tested that yet.
    It also binds some key EC_MS functions including limited (and un-tested)
    data importing (in __init__()), sync_metadata, add_selector, plot_vs_potential,
    and cutting (via cut_dataset or select_cycles) in cut().
    '''
    def __init__(self, file_name=None, folder=None, tag=None, data_type=None,
                 file_type=None, verbose=True):
        '''
        Establishes the dataset by loading self.data
        '''
        self.type = 'Dataset'
        if folder is not None: # then go to the folder and remember how to get back
            back = os.getcwd()
            os.chdir(folder)

        if type(file_name) is dict and 'data_cols' in file_name:
            #^ user can intiate the dataset with a data dictionary
            self.data = file_name
        elif type(file_name) in (list, tuple):
            # ^ user can intiate the dataset with a list of data files
            datas = []
            for file in file_name:
                data = get_data_from_file(file, verbose=verbose)
                datas += [data]
            self.data = synchronize(datas, verbose=verbose)
        elif file_name is not None:
            # ^ ...or just one data file
            self.data = get_data_from_file(file_name)
        elif folder is not None:
            # ^ ...or a bunch of files in a folder
            files = os.listdir() # note, we are already in the folder
            if tag is not None:
                files = [f for f in files if re.search('^' + tag, f)]
            if file_type is not None:
                files = [f for f in files if re.search(file_type + '$', f)]
            datas = []
            for file in files:
                data = get_data_from_file(file, verbose=verbose)
                datas += [data]
            self.data = synchronize(datas, verbose=verbose)
        else:
            print('Warning!!! Please specify file_name and/or folder.' +
                  ' Returning an empty dataset')

        if folder is not None: # time to go home.
            os.chdir(back)

        self.update_with_data()

    def update_with_data(self):
        if not hasattr(self, 'data') or 'data_cols' not in self.data:
            print('Warning!!! Empty dataset.')
            return
        for key, value in self.data.items():
            if key not in self.data['data_cols']:
                try:
                    setattr(self, key, value)
                except:
                    # not sure what the error is yet.
                    raise


    def __getattr__(self, attr):
        '''
        Makes it so that you can get items in self.data as if they were
        attributes to the Dataset.
        '''
        if attr == 't':
            return self.data[self.t_str]
        elif attr == 'v':
            return self.data[self.V_str]
        elif attr == 'j':
            return self.data[self.J_str]

        try:
            print('getting attribute ' + attr + ' from self.data') # debugging
            #a = b # debugging
            return self.data[attr]
        except KeyError:
            raise AttributeError('Dataset has no attribute ' + attr)

    def __getitem__(self, key):
        '''
        Makes it so that you can look up attributes to self as if they were
        items in a dictionary.
        Attributes pre-empt items in self.data.
        '''
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError('Dataset has no attribute ' + key +
                           ' and Dataset.data has no key ' + key
                           )

    def __setitem__(self, key, value):
        setattr(self, key, value)
        self.data[key] = value


    def __add__(self, dataset_2):
        new_data = synchronize([self.data, dataset_2.data])
        new_dataset = Dataset(new_data)
        return new_dataset



    def save(self, file_name, data_type=None):
        if data_type is None:
            data = self.data
        else:
            data = {}
            data_cols = set()
            for key, value in self.data.items():
                if key in self.data['data_cols']:
                    if self.data['col_types'][key] == data_type:
                        data_cols.add(key)
                        data[key] = value
                    else:
                        continue
                else:
                    data[key] = value
            data['data_cols'] = data_cols
            data['data_type'] = data_type
        with open(file_name, 'wb') as f:
            pickle.dump(data, f)

         # Binding existing functions. There is probably a much smarter way to do this...
    @wraps(sync_metadata)
    @with_update
    def sync_metadata(self, *args, **kwargs):
        print('args = ' + str(args)) # debugging. proves that args[0] is self.
        print('kwargs = ' + str(kwargs)) # debugging
        return sync_metadata(self.data, *args, **kwargs)

    @wraps(make_selector)
    @with_update
    def make_selector(self, *args, **kwargs):
        return make_selector(self.data, *args, **kwargs)

    @wraps(correct_ohmic_drop)
    @with_update
    def correct_ohmic_drop(self, *args, **kwargs):
        return correct_ohmic_drop(self.data, *args, **kwargs)

    @wraps(plot_experiment)
    def plot_experiment(self, *args, **kwargs):
        return plot_experiment(self.data, *args, **kwargs)

    @wraps(plot_vs_potential)
    def plot_vs_potential(self, *args, **kwargs):
        return plot_vs_potential(self.data, *args, **kwargs)

    # ... yes, there is! Just equate the function. If the getitem and getattr of
    # Dataset work as well as I hope, the function won't notice it's getting
    # the Dataset object as the first argument rather than the data dictionary.
    #sync_metadata = sync_metadata
    #make_selector = make_selector
    #correct_ohmic_drop = correct_ohmic_drop
    #plot_experiment = plot_experiment
    #plot_vs_potential = plot_vs_potential

    def normalize(self, *args, **kwargs):
        return self.sync_metadata(*args, **kwargs)
    def calibrate_EC(self, *args, **kwargs):
        return self.sync_metadata(*args, **kwargs)

    def cut(self, tspan=None, cycles=None, **kwargs):
        if tspan is not None:
            new_data = cut_dataset(self.data, tspan=tspan, **kwargs)
        else:
            for key in ['cycle number', 'selector', 'loop number',
                        'file number', 'cycle']:
                # should add self.sel_str, but that would require major changes
                if key in kwargs:
                    cycles = kwargs.pop(key)
                    new_data = select_cycles(self.data, cycles=cycles,
                                             cycle_str=key, **kwargs)
        new_dataset = Dataset(new_data)
        return new_dataset


class CyclicVoltammagram(Dataset):
    '''
    CyclicVoltammagram inherits from Dataset. It is easiest to initiate a CyclicVoltammagram
    by cutting a Dataset. The main addition is that it has a mandatory default
    selector called 'cycle', and indexing by this selects cycles.
    The default plotting function plot() is plot_vs_potential.
    It binds CV_difference(). It also has a few brand new functions including
    redefine_cycle() which lets you say where CVs start,
    plot_all() which plots cv's with a cmap, average() with averages cycles.
    '''

    def __init__(self, *args, **kwargs):
        self.type = 'Cyclic Voltammagram'
        if 'dataset' in kwargs:
            dataset = kwargs.pop('dataset')
        else:
            dataset = args[0] # 'tuple' object has no attribute 'pop' :(
            if len(args)>0:
                args = args[1:]
            else:
                args = []

        if type(dataset) is dict:
            dataset = Dataset(dataset)
        elif hasattr(dataset, 'type'):
            if dataset.type == 'Dataset':
                if len(args)>0 or len(kwargs)>0:
                    kwargs.update(verbose=False)
                    dataset = dataset.cut(*args, **kwargs)
                data = dataset.data
            else:
                print('WARNING!!! CyclicVoltammagram.__init__ doesn\'t know ' +
                      'what ' + str(dataset) + ' is!')
                data = {}
            self.data = data
            self.update_with_data()
            self.redefine_cycle()
        else:
            dataset = Dataset(dataset, *args, **kwargs)
            self.__init__(dataset)


    def __getitem__(self, key):
        '''
        Makes it so that you can look up attributes to self as if they were
        items in a dictionary.
        Attributes pre-empt items in self.data.
        '''
        if type(key) is slice:
            start, stop, step = key.start, key.stop, key.step
            if step is None:
                step = 1
            key = list(range(start, stop, step))
        if type(key) in [int, list]:
            return CyclicVoltammagram(self.cut(cycle=key, t_zero='start',
                                               verbose=False))
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError('Dataset has no attribute ' + key +
                           'and Dataset.data has no key ' + key
                           )

    def __len__(self):
        return len(set(self.cycle))

    def redefine_cycle(self, V=None, redox=None):
        if V is None:
            try:
                selector = self[self['sel_str']]
            except KeyError:
                sel_str = self.make_selector()
                selector = self[sel_str]
            cycle = selector - min(selector)

        else:
            cycle = np.zeros(self.t.shape)
            c = 0
            n = 0
            N = len(self.t)
            v = self.v
            if redox in [0, -1, 'red', 'reduction']:
                # easiest way to reverse directions is to use the same > < operators
                # but negate the arguments
                V = -V
                v = -v
            while n<N:
                mask_behind = v[n:]<V
                if not True in mask_behind:
                    break
                else:
                    n += np.argmax(mask_behind) + 5 # have to be below V for 5 datapoints
                print('point number on way up: ' + str(n)) # debugging

                mask_in_front = v[n:]>V
                if not True in mask_in_front:
                    break
                else:
                    n += np.argmax(mask_in_front)
                c += 1 # and then when it crosses to above V again, we register a cyclce!
                cycle[n:] = c # and subsequent points increase in cycle number
                n += + 5 # have to be above V for 5 datapoints
                print('point number on way down: ' + str(n)) # debugging

        self['cycle'] = cycle
        self.data['sel_str'] = 'cycle'
        self.data_cols.add('cycle')


    @wraps(plot_vs_potential)
    def plot(self, *args, **kwargs):
        return self.plot_vs_potential(*args, **kwargs)


    @wraps(CV_difference)
    def get_difference(self, *args, **kwargs):
        return CV_difference(self.data, *args, **kwargs)


    def plot_all(self, ax='new', colorscale='spectral', **kwargs):
        cmap = colormap.get_cmap(colorscale)

        if ax == 'new':
            fig, ax = plt.subplots()

        C = len(self)
        for c in range(C):
            color = cmap(c/C)
            self[c].plot(ax=ax, color=color, **kwargs)

        ax.set_xlabel(self.V_str)
        ax.set_ylabel(self.J_str)

        return ax

    def average(self):

        # should in principle do this for much more

        lists = {}
        for col in self.data_cols:
            lists[col] = []
        for c in range(len(self)):
            cv = self[c]
            for col in self.data_cols:
                lists[col] += [cv[col]]
        ts = lists['time/s']
        N = min([len(t) for t in ts])

        average_cv = self[0] # inherit all the metadata from self[0]
        for col in self.data_cols:
            x_stack = np.stack([x[:N] for x in lists[col]])
            x_avg = np.mean(x_stack, axis=0)

            average_cv.data[col] = x_avg

        average_cv.update_with_data()
        return average_cv





