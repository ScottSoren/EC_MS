# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 19:08:46 2020

@author: scott
"""
import os, pickle
import numpy as np
from functools import wraps
from matplotlib import pyplot as plt
from .dataset import Dataset
from .spectra import Spectrum, Spectra, spectra_from_data

"""
The main Zilien importing is at present taken care of by .Data_Importing/load_from_file
and the chaotic multi-format parser that it calls.
I think a better way would be to have a module for each data type, inhereting from Dataset
and with its own parsers, which may use some shared tools in a shared module.
In general, the structure of EC_MS needs serious reworking!
"""


class Zilien_Dataset(Dataset):
    # @wraps(Dataset.__init__)
    def __init__(self, *args, **kwargs):
        if "data_type" not in kwargs:
            kwargs["data_type"] = "SI"
        print(kwargs)
        super().__init__(*args, **kwargs)
        self.get_spectra()

    def get_spectra(self):
        if "spectra_data" in self.data:
            self.spectra = spectra_from_data(self.data["spectra_data"])
        else:
            try:
                spectra_folder = (
                    "".join([s + " " for s in self.file.split(" ")[2:]]).split(".")[0]
                    + " mass scans"
                )
                spectra_path = os.path.join(self.folder, spectra_folder)
                self.spectra_folder, self.spectra_path = spectra_folder, spectra_path
                self.spectra = read_zilien_spectra(spectra_path)
            except FileNotFoundError:
                print("Warning!!! No spectra found! consider using normal Dataset")
            # self.spectrums = self.spectra.spectrums

    def __getitem__(self, key):
        if type(key) is int:
            return self.spectra[key]

    def save(self, file_name):
        spectra_data = {"x": self[0].x, "spectra": self.spectra.spectra}
        self.data["spectra_data"] = spectra_data
        with open(file_name, "wb") as f:
            pickle.dump(self.data, f)


def read_zilien_spectrum(file_path, delim="\t"):

    with open(file_path, "r") as f:
        lines = f.readlines()

    data = {
        "file": file_path,
        "header": "",
    }
    N_col_head = len(
        lines
    )  # this will decrease when the loop knows when the column header line is comming

    nondata_cols = []  # this will store abstime, so I don't have to parse.
    for n, line in enumerate(lines):
        l = line.strip()
        if n < N_col_head:
            if len(l) == 0:
                N_col_head = n + 1
            # print(dataset['header']) # debugging
            # if n< 10: print(line)
            # data['header'] = data['header'] + line # If I use .join instead, it gives a memory error, I don't understand why.
        elif n == N_col_head:
            data_cols = l.split(delim)
            for col in data_cols:
                data[col] = np.array([])
            # data['header'] = data['header'] + line # If I use .join instead, it gives a memory error, I don't understand why.
        elif n > N_col_head:
            for col, val in zip(data_cols, l.split(delim)):
                if col in nondata_cols:
                    data[col] += [val]
                    continue
                try:
                    x = eval(val)
                except SyntaxError:
                    print(
                        "removing "
                        + col
                        + " from data_cols due to value "
                        + val
                        + " on line "
                        + str(n)
                    )
                    data[col] = list(data[col])
                    data[col] += [val]
                    nondata_cols += [col]
                else:
                    data[col] = np.append(data[col], x)

    data["data_cols"] = set(data_cols)

    return data


def read_zilien_spectrums(folder, delim="\t"):
    """
    """
    lslist = os.listdir(folder)
    spectra = []
    ts = []
    for f in lslist:
        try:
            time_str = f.split("started at measurement time")[1]
        except IndexError:
            print(f + " does not seem to be a Zilien spectrum with timestamp")
        else:
            time_str = time_str.split(".tsv")[0].strip()
        t = float(time_str)
        data = read_zilien_spectrum(folder + os.sep + f, delim=delim)
        # return data # debugging
        x = data["Mass  [AMU]"]
        y = data["Current [A]"]
        spectrum = Spectrum(x=x, y=y, t=t)
        ts += [t]
        spectra += [spectrum]
    I_sort = np.argsort(ts)
    ts = [ts[I] for I in I_sort]
    spectrums = [
        spectra[I] for I in I_sort
    ]  # can't directly write spectra[I_sort] since it's not an np array
    return spectrums


def read_zilien_spectra(folder, delim="\t"):
    """
    """
    spectrums = read_zilien_spectrums(folder, delim="\t")
    # return spectrums # debugging
    return Spectra(folder=folder, spectrums=spectrums)
