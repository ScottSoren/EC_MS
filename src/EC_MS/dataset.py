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
from .Combining import synchronize, cut_dataset, sort_time, get_timecol
from .Plotting import plot_experiment, plot_vs_potential, plot_flux, plot_signal
from .EC import correct_ohmic_drop, CV_difference
from .Quantification import get_current, get_signal, get_potential
from .Calibration import calibration_curve, point_calibration, chip_calibration


def get_data_from_file(
    file_name, data_type=None, verbose=True
):  # assumes you're already in the folder
    if type(file_name) is dict:
        return file_name  # so that the dataset can be initiated with data already in a dictionary
    if re.search(".pkl$", file_name):
        with open(file_name, "rb") as f:
            return pickle.load(f)
    elif data_type is not None:
        return load_from_file(file_name, data_type=data_type, verbose=verbose)
    elif re.search(".mpt$", file_name):
        return load_from_file(file_name, data_type="EC", verbose=verbose)
    else:
        print(
            "WARNING: loading files of the type "
            + file_name
            + " is not yet implemented in Dataset.__init__() !!!"
            + " Try specifying a data_type."
        )


def with_update(method):
    @wraps(method)
    def method_with_update(self, *args, **kwargs):
        ret = method(self, *args, **kwargs)
        self.update_with_data()
        return ret

    return method_with_update


metadata_items = [
    "data_type",
    "mass_bgs",
]


class Dataset:
    """
    This class implements the dataset. Its design is to be back-compatible
    with the dataset dictionaries that were the main object in function-
    centric EC_MS programming.

    Dataset just serves as a wrapper around dataset dictionaries to
    make the package seem object-oriented. It has __getitem__ and
    __getattr__ methods that make (key,value) pairs and attributes somewhat
    interchangable. It should be back-compatable, but I haven't really tested that yet.
    It also binds some key EC_MS functions including various options for data importing in
    data importing in __init__(); scutting (via cut_dataset or select_cycles) in cut();
    plus several plotting (plot_experiment, plot_vs_potential) and getting (get_signal, get_current) functions
    """

    def __init__(
        self,
        file_name=None,
        folder=None,
        tag=None,
        data_type=None,
        file_type=None,
        verbose=True,
    ):
        """
        Establishes the dataset by loading self.data
        """
        self.type = "Dataset"
        self.verbose = verbose
        if folder is not None:  # then go to the folder and remember how to get back
            back = os.getcwd()
            os.chdir(folder)

        if type(file_name) is dict and "data_cols" in file_name:
            # ^ user can intiate the dataset with a data dictionary
            self.data = file_name
        elif type(file_name) in (list, tuple):
            # ^ user can intiate the dataset with a list of data files
            datas = []
            for file in file_name:
                data = get_data_from_file(file, verbose=verbose, data_type=data_type)
                datas += [data]
            self.data = synchronize(datas, verbose=verbose)
        elif file_name is not None:
            # ^ ...or just one data file
            self.data = get_data_from_file(file_name, data_type=data_type)
        elif folder is not None:
            # ^ ...or a bunch of files in a folder
            print("Importing from a folder!!!")
            files = os.listdir()  # note, we are already in the folder
            if tag is not None:
                files = [f for f in files if re.search("^" + tag, f)]
            if file_type is not None:
                files = [f for f in files if re.search(file_type + "$", f)]
            # print(files) # debugging
            datas = []
            for file in files:
                data = get_data_from_file(file, verbose=verbose)
                datas += [data]
            self.data = synchronize(datas, verbose=verbose)
            sort_time(self.data)
        else:
            print(
                "Warning!!! Please specify file_name and/or folder."
                + " Returning an empty dataset"
            )

        if folder is not None:  # time to go home.
            os.chdir(back)
        if data_type is not None:
            self.data["data_type"] = data_type
        self.update_with_data()

    def update_with_data(self):
        if not hasattr(self, "data") or "data_cols" not in self.data:
            print("Warning!!! Empty dataset.")
            return
        for key, value in self.data.items():
            if key not in self.data["data_cols"]:
                try:
                    setattr(self, key, value)
                except:
                    # not sure what the error is yet.
                    raise

    def __getattr__(self, attr):
        """
        Makes it so that you can get items in self.data as if they were
        attributes to the Dataset.
        """
        if attr == "t":
            return self.data[self.t_str]
        elif attr == "v":
            return self.data[self.V_str]
        elif attr == "j":
            return self.data[self.J_str]
        elif attr in ["data", "verbose"]:
            # print('hm...' + attr) # debugging
            raise AttributeError("Dataset has no attribute " + attr)
        # print('getting attribute ' + attr + ' from self.data') # debugging
        try:
            # a = b # debugging
            return self.data[attr]
        except KeyError:
            raise AttributeError("Dataset has no attribute " + attr)

    def __getitem__(self, key):
        """
        Makes it so that you can look up attributes to self as if they were
        items in a dictionary.
        Attributes pre-empt items in self.data.
        """
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(
                "Dataset has no attribute "
                + key
                + " and Dataset.data has no key "
                + key
            )

    def __setitem__(self, key, value):
        setattr(self, key, value)
        self.data[key] = value

    def __add__(self, dataset_2):
        new_data = synchronize([self.data, dataset_2.data])
        new_dataset = Dataset(new_data)
        return new_dataset

    def add_data_col(self, col, value, col_type=None):
        self.data[col] = value
        self.data["data_cols"].add(col)
        if col_type is not None:
            if "col_types" not in self.data:
                self.data["col_types"] = {}
            self.data["col_types"][col] = col_type

    def append_to_data_col(self, col, value, col_type=None):
        """
        Handy thing, adds data to a col if it exist, otherwise creates the
        col and puts it in self.data_cols
        """
        if col in self.data_cols:
            self.data[col] = np.append(self.data[col], value)
        else:
            self.add_data_col(col, value, col_type)

    def save(self, file_name, data_type=None):
        if data_type is None:
            data = self.data
        else:
            data = {}
            data_cols = set()
            for key, value in self.data.items():
                if key in self.data["data_cols"]:
                    if self.data["col_types"][key] == data_type:
                        data_cols.add(key)
                        data[key] = value
                    else:
                        continue
                else:
                    data[key] = value
            data["data_cols"] = data_cols
            data["data_type"] = data_type
        with open(file_name, "wb") as f:
            pickle.dump(data, f)

        # Binding existing functions. There is probably a much smarter way to do this...

    @wraps(sync_metadata)
    @with_update
    def sync_metadata(self, *args, **kwargs):
        # print('args = ' + str(args)) # debugging. proves that args[0] is self.
        # print('kwargs = ' + str(kwargs)) # debugging
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

    @wraps(plot_flux)
    def plot_flux(self, *args, **kwargs):
        return plot_flux(self.data, *args, **kwargs)

    @wraps(plot_signal)
    def plot_signal(self, *args, **kwargs):
        return plot_signal(self.data, *args, **kwargs)

    @wraps(plot_vs_potential)
    def plot_vs_potential(self, *args, **kwargs):
        return plot_vs_potential(self.data, *args, **kwargs)

    @wraps(get_current)
    def get_current(self, *args, **kwargs):
        return get_current(self.data, *args, **kwargs)

    @wraps(get_signal)
    def get_signal(self, *args, **kwargs):
        return get_signal(self.data, *args, **kwargs)

    @wraps(get_potential)
    def get_potential(self, *args, **kwargs):
        return get_potential(self.data, *args, **kwargs)

    @wraps(get_timecol)
    def get_timecol(self, col, **kwargs):
        return get_timecol(col, dataset=self.data, **kwargs)

    @wraps(calibration_curve)
    def calibration_curve(self, **kwargs):
        return calibration_curve(self.data, **kwargs)

    @wraps(point_calibration)
    def point_calibration(self, **kwargs):
        return point_calibration(self.data, **kwargs)

    @wraps(chip_calibration)
    def chip_calibration(self, **kwargs):
        return chip_calibration(self.data, **kwargs)

    def get_flux(self, m, *args, **kwargs):
        try:
            return m.get_flux(self.data, *args, **kwargs)
        except AttributeError:
            print(
                "WARNING!!! first argument to dataset.get_flux must be an object of class EC_MS.Molecule."
            )
            raise TypeError

    # ... yes, there is! Just equate the function. If the getitem and getattr of
    # Dataset work as well as I hope, the function won't notice it's getting
    # the Dataset object as the first argument rather than the data dictionary.
    # sync_metadata = sync_metadata
    # make_selector = make_selector
    # correct_ohmic_drop = correct_ohmic_drop
    # plot_experiment = plot_experiment
    # plot_vs_potential = plot_vs_potential

    def normalize(self, *args, **kwargs):
        return self.sync_metadata(*args, **kwargs)

    def calibrate_EC(self, *args, **kwargs):
        return self.sync_metadata(*args, **kwargs)

    def set_background(self, t_bg=None, masses=None, mols=None, cols=None):
        """
        TODO: if given mols, it sets a background in each of of the given molecule
        objects using mol.get_bg #ToDo, that should be mol.set_bg instead

        if given masses, it calculates the background of each
        and stores it in a dictionary, and SUBTRACTS IT FROM THE DATA!

        TODO: don't subtract it from the data, but have get_signal read it.
        """
        self.reset()  # to avoid losing the ability to restore the original
        # by subtracting a new background from background-subtracted data

        if masses is None and mols is None and cols is None and t_bg is not None:
            masses = "all"
        if masses == "all":
            masses = [
                col[:-2]
                for col in self.data_cols
                if (col[0] == "M" and col[-2:] == "-y")
            ]
        print("masses = " + str(masses))  # debugging

        if hasattr(self, "mass_bgs"):
            mass_bgs = self.mass_bgs
        else:
            mass_bgs = {}

        if masses is not None:
            for mass in masses:
                if t_bg is not None:
                    x, y = self.get_signal(mass=mass, tspan=t_bg, unit="A")
                    y_bg = np.mean(y)
                else:
                    x, y = self.get_signal(mass=mass, unit="A")
                    y_bg = min(y)
                # print('subtracting background for mass ' + mass + '!!!')
                mass_bgs[mass] = y_bg
                self.data[mass + "-y"] -= y_bg
        self.mass_bgs = mass_bgs
        return mass_bgs

    def reset(self):
        """
        so far only implemented for masses.
        """
        if hasattr(self, "mass_bgs"):
            for mass, y_bg in self.mass_bgs.items():
                print("adding background back onto " + mass)  # debugging
                self.data[mass + "-y"] += y_bg

    def cut(self, tspan=None, cycles=None, verbose=True, **kwargs):
        if tspan is not None:
            new_data = cut_dataset(self.data, tspan=tspan, **kwargs)
        else:
            for key in [
                "cycle number",
                "selector",
                "loop number",
                "file number",
                "cycle",
                "sweep",
            ]:
                # should add self.sel_str, but that would require major changes
                if key in kwargs:
                    cycles = kwargs.pop(key)
                    new_data = select_cycles(
                        self.data,
                        cycles=cycles,
                        cycle_str=key,
                        verbose=verbose,
                        **kwargs
                    )
        new_dataset = Dataset(new_data)
        for attr in metadata_items:
            if hasattr(self, attr):
                setattr(new_dataset, attr, getattr(self, attr))
        return new_dataset


class CyclicVoltammagram(Dataset):
    """
    CyclicVoltammagram inherits from Dataset. It is easiest to initiate a CyclicVoltammagram
    by cutting a Dataset. The main addition is that it has a mandatory default
    selector called 'cycle', and indexing by this selects cycles.
    The default plotting function plot() is plot_vs_potential.
    It binds CV_difference(). It also has a few brand new functions including
    redefine_cycle() which lets you say where CVs start,
    plot_all() which plots cv's with a cmap, average() with averages cycles.
    """

    def __init__(self, *args, verbose=True, **kwargs):
        self.type = "Cyclic Voltammagram"
        self.verbose = verbose
        if "dataset" in kwargs:
            dataset = kwargs.pop("dataset")
        else:
            dataset = args[0]  # 'tuple' object has no attribute 'pop' :(
            if len(args) > 0:
                args = args[1:]
            else:
                args = []

        if type(dataset) is dict:
            dataset = Dataset(dataset)

        if hasattr(dataset, "type"):
            if dataset.type in ["Dataset", "Cyclic Voltammagram"]:
                if len(args) > 0 or len(kwargs) > 0:
                    kwargs.update(verbose=False)
                    dataset = dataset.cut(*args, **kwargs)
                data = dataset.data
                for attr in metadata_items:
                    if hasattr(dataset, attr):
                        setattr(self, attr, getattr(dataset, attr))
            else:
                print(
                    "WARNING!!! CyclicVoltammagram.__init__ doesn't know "
                    + "what "
                    + str(dataset)
                    + " is!"
                )
                data = {}
            self.data = data
            self.update_with_data()
            self.redefine_cycle()
        else:
            dataset = Dataset(dataset, *args, **kwargs)
            self.__init__(dataset)

    def __getitem__(self, key):
        """
        Makes it so that you can look up attributes to self as if they were
        items in a dictionary.
        Attributes pre-empt items in self.data.
        """
        if type(key) is slice:
            start, stop, step = key.start, key.stop, key.step
            if step is None:
                step = 1
            key = list(range(start, stop, step))
        if type(key) in [int, list]:
            return CyclicVoltammagram(
                self.cut(cycle=key, t_zero="start", verbose=False), verbose=False
            )
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(
                "Dataset has no attribute "
                + key
                + " and Dataset.data has no key "
                + key
            )

    def __len__(self):
        return len(set(self.cycle))

    def redefine_cycle(self, V=None, redox=None):
        """
        Changes self.data['cycle'] to count each time the calibrated potential
        passes V in the direction specified by redox (1 for anodic, 0 for cathodic)
        """
        if V is None:
            try:
                selector = self[self["sel_str"]]
            except KeyError:
                sel_str = self.make_selector()
                # print(self.data.keys()) # debugging
                selector = self[sel_str]
            cycle = selector - min(selector)

        else:
            cycle = np.zeros(self.t.shape)
            c = 0
            n = 0
            N = len(self.t)
            v = self.v
            if redox in [0, -1, "red", "reduction"]:
                # easiest way to reverse directions is to use the same > < operators
                # but negate the arguments
                V = -V
                v = -v
            while n < N:
                mask_behind = v[n:] < V
                if not True in mask_behind:
                    break
                else:
                    n += (
                        np.argmax(mask_behind) + 5
                    )  # have to be below V for 5 datapoints
                # print('point number on way up: ' + str(n)) # debugging

                mask_in_front = v[n:] > V
                if not True in mask_in_front:
                    break
                else:
                    n += np.argmax(mask_in_front)
                c += 1  # and then when it crosses to above V again, we register a cyclce!
                cycle[n:] = c  # and subsequent points increase in cycle number
                n += +5  # have to be above V for 5 datapoints
                # print('point number on way down: ' + str(n)) # debugging

        self.add_data_col("cycle", cycle, col_type="EC")
        self.data["sel_str"] = "cycle"

    def get_sweeps(self, min_sweep_points=10, scan_rate_cutoff=1):
        try:
            return self.data["sweep_types"]
        except KeyError:
            return self.make_sweeps(
                min_sweep_points=min_sweep_points, scan_rate_cutoff=scan_rate_cutoff
            )

    def make_sweeps(self, min_sweep_points=10, scan_rate_cutoff=1):
        """
        figures out when anodic (sweep_type=1) and cathodic (sweep_type=0) sweeps
        and potential holds (sweep_type=None) are in the data based on the
        scan rate.
        min_sweep_points is the resolution in EC points.
        scan_rate_cutoff is the minimum absolute scan rate in mV needed to be considered
        an anodic or cathodic sweep
        """
        print("\n\nfunction CyclicVoltammagram.make_sweeps at your service!\n")
        sweep_types = {}  # 1 for oxidation, -1 for reduction, None for staying still
        sweep_index_to_sweep_type = {
            0: 0,
            1: 1,
            2: None,
        }  # but we'll use a list in the grind
        sweep = np.zeros(self.t.shape)

        scan_rate = self.get_scan_rate(min_sweep_points=min_sweep_points)

        cat_mask = scan_rate < -scan_rate_cutoff
        an_mask = scan_rate > scan_rate_cutoff
        hold_mask = abs(scan_rate) < scan_rate_cutoff

        the_masks = [an_mask, cat_mask, hold_mask]
        for mask in the_masks:
            mask[
                -2
            ] = False  # because np.argmin(mask)=0 if mask is True all the time, giving problems
            mask[
                -1
            ] = True  # because np.argmax(mask)=0 if mask is False all the time, giving problems

        # print('the_masks:\n' + str(the_masks)) # debugging

        N = len(self.t)
        i_start = 0
        i_finish = 0
        n_sweep = 0

        the_next_starts = [np.argmax(mask) for mask in the_masks]
        sweep_index = np.argmin(the_next_starts)

        while i_start < N - 1:
            # print('\n\n') # debugging
            # print('the next starts = ' + str(the_next_starts)) # debugging
            I_out = np.argmin(the_masks[sweep_index][i_finish:])
            # print(the_masks[sweep_index][i_finish:]) # debugging
            i_start = i_finish + I_out + min_sweep_points
            # can't start a new sweep until you've been out of the current sweep for at least min_sweep_points

            try:
                I_in_again = np.argmax(the_masks[sweep_index][i_start:])
            except ValueError:
                the_next_starts[sweep_index] = N
            else:
                # print(the_masks[sweep_index][i_start:]) # debugging
                # ^ check how long until the next sweep of that type starts
                the_next_starts[sweep_index] = i_start + I_in_again
                # ^ and add it.

            next_sweep_index = np.argmin(the_next_starts)
            i_finish = the_next_starts[next_sweep_index]

            # print('I_out = ' + str(I_out) + ', I_in_again = ' + str(I_in_again)) # debugging
            # print('i_start = ' + str(i_start) + ', i_finish = ' + str(i_finish)) # debugging
            # print('sweep index = ' + str(sweep_index) + ', the_masks[sweep_index] = ' + str(the_masks[sweep_index]))
            # if n_sweep > 10: break # debugging

            if not next_sweep_index == sweep_index:
                sweep_index = next_sweep_index
                sweep[i_finish:] += 1
                sweep_types[n_sweep] = sweep_index_to_sweep_type[sweep_index]
                n_sweep += 1

        self.add_data_col("sweep", sweep, "EC")
        self.data["sweep_types"] = sweep_types
        print("\nfunction CyclicVoltammagram.make_sweeps finished!\n\n")
        return self.data["sweep_types"]

    def get_scan_rate(self, min_sweep_points=10, tspan=None, cycle=None):
        """
        returns scan rate in mV/s. If a tspan or cycle is given, it returns
        the average absolute scan rate for that time interval or cycle.
        Otherwise it returns a vector.
        """
        try:
            scan_rate = self.data["scan_rate"]
        except KeyError:
            scan_rate = self.make_scan_rate(min_sweep_points=min_sweep_points)

        if cycle is not None:
            return np.mean(np.abs(self[cycle]).scan_rate)
        if tspan is not None:
            t, scan_rate = self.t, self.scan_rate
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            return np.mean(scan_rate[mask])
        return self.scan_rate

    def make_scan_rate(self, min_sweep_points=10):
        """
        calculates scan rate in mV/s - negative for cathodic, positive for anodic.
        min_sweep_points is a type of resolution in EC points.
        """
        print("\n\nfunction CyclicVoltammagram.make_scan_rate at your service!\n")
        v = self.v
        t = self.t

        # the scan rate is dV/dt. This is a numerical calculation of dV/dt:
        v_behind = np.append(np.tile(v[0], min_sweep_points), v[:-min_sweep_points])
        v_ahead = np.append(v[min_sweep_points:], np.tile(v[-1], min_sweep_points))

        t_behind = np.append(np.tile(t[0], min_sweep_points), t[:-min_sweep_points])
        t_ahead = np.append(t[min_sweep_points:], np.tile(t[-1], min_sweep_points))

        scan_rate_middle = (v_ahead - v_behind) / (t_ahead - t_behind) * 1e3
        # ^ this is "softened" at the anodic and cathodic turns.

        # We can "sharpen" it by selectively looking ahead and behind:
        scan_rate_behind = (v - v_behind) / (t - t_behind) * 1e3
        scan_rate_ahead = (v_ahead - v) / (t_ahead - t) * 1e3

        # but this gives problems right at the beginning, so set those to zeros
        scan_rate_behind[:min_sweep_points] = np.zeros(min_sweep_points)
        scan_rate_ahead[-min_sweep_points:] = np.zeros(min_sweep_points)

        # now sharpen the scan rate!
        scan_rate = scan_rate_middle
        mask_use_ahead = np.logical_and(
            np.abs(scan_rate_ahead) > np.abs(scan_rate),
            np.abs(scan_rate_ahead) > np.abs(scan_rate_behind),
        )
        scan_rate[mask_use_ahead] = scan_rate_ahead[mask_use_ahead]
        mask_use_behind = np.logical_and(
            np.abs(scan_rate_behind) > np.abs(scan_rate),
            np.abs(scan_rate_behind) > np.abs(scan_rate_ahead),
        )
        scan_rate[mask_use_behind] = scan_rate_behind[mask_use_behind]

        if False:  # plot it for debugging
            fig, ax = plt.subplots()
            ax.plot(t, scan_rate, "k")  # t_ahead - t_behind,

        self.add_data_col("scan_rate", scan_rate, "EC")
        print("\nfunction CyclicVoltammagram.make_scan_rate finished!\n\n")
        return self.scan_rate

    @wraps(plot_vs_potential)
    def plot(self, *args, **kwargs):
        return self.plot_vs_potential(*args, **kwargs)

    @wraps(CV_difference)
    def get_difference(self, *args, **kwargs):
        return CV_difference(self.data, *args, **kwargs)

    def subtract(self, cv_2, min_sweep_points=10, scan_rate_cutoff=1):
        """
        takes the datasets sweep by sweep, interpolates, and
        """
        print("\n\nfunction CyclicVoltammagram.subtract at your service!\n")
        # best to force a remake of the sweep numbers.
        sweeps_1 = self.make_sweeps(
            min_sweep_points=min_sweep_points, scan_rate_cutoff=scan_rate_cutoff
        )
        sweeps_2 = cv_2.make_sweeps(
            min_sweep_points=min_sweep_points, scan_rate_cutoff=scan_rate_cutoff
        )

        s1s = list(sweeps_1.keys())
        for s in s1s:
            if sweeps_1[s] is None:
                s1s.pop(s)
        s1s.sort()
        s2s = list(sweeps_2.keys())
        for s in s2s:
            if sweeps_2[s] is None:
                s2s.pop(s)
        s2s.sort()

        try:
            t_str = self.t_str
        except AttributeError:
            t_str = "time/s"
        try:
            V_str = self.V_str
        except AttributeError:
            V_str = "Ewe/V"
        try:
            J_str = self.J_str
        except AttributeError:
            J_str = "I/mA"

        diff = Dataset(
            {"data_cols": set(), "V_str": V_str, "t_str": t_str, "J_str": J_str,}
        )
        debugging = False
        if debugging:
            fig, ax = plt.subplots()

        for s1, s2 in zip(s1s, s2s):
            print(
                "interpolating sweep "
                + str(s2)
                + " of cv_2 (redox = "
                + str(sweeps_2[s2])
                + ") onto sweep "
                + str(s1)
                + " of self (redox = "
                + str(sweeps_1[s1])
                + ")."
            )
            if not sweeps_2[s2] == sweeps_1[s1]:
                print(
                    "WARNING!!! subtracting sweeps of different directions. Results may be meaningless."
                )

            data1 = self.cut(sweep=s1).data
            data2 = cv_2.cut(sweep=s2).data

            t1, v = data1[t_str], data1[V_str]
            t2_i, v2 = data2[t_str], data2[V_str]

            if sweeps_1[s1] == 1:  # anodic scan, interpolate normally
                t2 = np.interp(v, v2, t2_i)
            elif sweeps_1[s1] == 0:  # cathodic scan, have to reverse the sine of v.
                t2 = np.interp(-v, -v2, t2_i)

            for col in data1["data_cols"]:

                if col in [t_str, V_str]:
                    diff.append_to_data_col(col, data1[col], col_type="EC")
                    # this will add it if it's not already there.
                elif col in data2["data_cols"]:
                    tcol = self.get_timecol(col)
                    if col == tcol:
                        continue
                    x1, y1_i = data1[tcol], data1[col]
                    x2, y2_i = data2[tcol], data2[col]

                    y1 = np.interp(t1, x1, y1_i)
                    y2 = np.interp(t2, x2, y2_i)

                    y_diff = y1 - y2

                    if debugging:
                        if "J" in col:
                            ax.plot(t1, y1, "k")
                            ax.plot(t1, y2, "r")
                            ax.plot(t1, y_diff, "g")

                    diff.append_to_data_col(
                        col, y_diff, col_type="EC"
                    )  # make it all bloody EC

        diff["tspan"] = [diff[t_str][0], diff[t_str][-1]]
        diff = CyclicVoltammagram(diff)
        try:
            diff.data["data_type"] = self.data["data_type"]
        except KeyError:
            print("Warning!!! no self.data['data_type']")

        print("\nfunction CyclicVoltammagram.subtract finished!\n\n")
        return diff

    def plot_all(self, ax="new", colorscale="spectral", **kwargs):
        cmap = colormap.get_cmap(colorscale)

        if ax == "new":
            fig, ax = plt.subplots()

        C = len(self)
        for c in range(C):
            color = cmap(c / C)
            self[c].plot(ax=ax, color=color, **kwargs)

        ax.set_xlabel(self.V_str)
        ax.set_ylabel(self.J_str)

        return ax

    def average(self):

        # should in principle do this for much more

        lists = {}
        for col in self.data_cols:
            lists[col] = []
        # print('test') # debugging
        for c in range(len(self)):
            cv = self[c]
            for col in self.data_cols:
                lists[col] += [cv.data[col]]
        ts = lists["time/s"]
        N = min([len(t) for t in ts])

        average_cv = self[0]  # inherit all the metadata from self[0]
        for col in self.data_cols:
            x_stack = np.stack([x[:N] for x in lists[col]])
            x_avg = np.mean(x_stack, axis=0)

            average_cv.data[col] = x_avg

        average_cv.update_with_data()
        return average_cv

    def cut(self, *args, **kwargs):
        dataset = super(CyclicVoltammagram, self).cut(*args, **kwargs)
        return CyclicVoltammagram(dataset)
