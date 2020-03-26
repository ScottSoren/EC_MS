# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:55:31 2020

@author: scott
"""
import os, re, pickle
import numpy as np
from matplotlib import pyplot as plt

from .PVMassSpec import read_PVMS

float_match = "[-]?\d+[\.]?\d*(e[-+]?\d+)?"  # matches floats like '-3.5e4' or '7' or '245.13' or '1e-15'


class Spectrum:
    def __init__(
        self,
        file_path=None,
        x=None,
        y=None,
        t=None,
        tstamp=None,
        index=0,
        data_type="PVMS",
    ):
        print("Initiating Spectrum!")  # debugging
        if x is None and y is None:
            spectra = Spectra(file_path, data_type=data_type)
            spectrum = spectra[index]
            x, y, t, tstamp = spectrum.x, spectrum.y, spectrum.t, spectrum.tstamp
        print(file_path)  # debugging
        self.x = x
        self.y = y
        self.t = t
        self.tstamp = tstamp
        self.file_path = file_path
        self.bg = 0
        if file_path is not None:
            self.folder, self.file = os.path.split(file_path)

    def get_signal(self, Mspan=None):
        x, y = self.x, self.y
        if Mspan is not None:
            mask = np.logical_and(Mspan[0] < x, x < Mspan[-1])
            x, y = x[mask], y[mask]
        return x, y

    def plot(self, Mspan=None, ax="new", **kwargs):
        """
        plots the spectrum. kwargs are fed to matplotlib.pyplot.plot
        """
        if ax == "new":
            fig, ax = plt.subplots()
            ax.set_xlabel("m/z")
            ax.set_ylabel("signal")
        x, y = self.get_signal(Mspan=Mspan)
        ax.plot(x, y, **kwargs)
        return ax

    def reset(self):
        """
        undoes self.set_background
        """
        self.y = self.y + self.bg
        self.bg = 0

    def set_background(self, M_bg=None, Mspan=None):
        """
        sets self.bg to the average value within Mspan.
        If Mspan is none, sets self.bg to min(self.y)
        Either way, saves self.bg and subtracts it from self.y
        self.reset() undoes this.
        This function calls self.reset() at the beginning to not compound
        background sets
        """
        self.reset()
        if M_bg is None and Mspan is not None:  # I can see this being a mistake
            M_bg = Mspan
        if M_bg is not None:
            x_bg, y_bg = self.get_signal(Mspan=M_bg)
            bg = np.mean(y_bg)
        else:
            bg = min(self.y)
        self.y = self.y - bg
        self.M_bg = M_bg
        self.bg = bg
        return bg

    def get_max(self, Mspan=None, ax=None):
        x, y = self.get_signal(Mspan=Mspan)
        return max(y)

    def get_integral(self, Mspan=None, M_bg=None, ax=None, **kwargs):
        x, y = self.get_signal(Mspan=Mspan)
        if M_bg is not None:
            x_bg, y_bg = self.get_signal(Mspan=M_bg)
            bg = np.mean(y_bg)
        else:
            bg = 0
        integral = np.trapz(y - bg, x)
        if ax is not None:
            if ax == "new":
                ax = self.plot(Mspan=Mspan)
            ax.fill_between(x, y, bg * np.ones(x.shape), **kwargs)
        return integral


def spectrums_from_data(data):
    x = data["x"]
    spectra = data["spectra"]
    spectrums = []
    for i, y in enumerate(spectra):
        spectrum = Spectrum(x=x, y=y)
        spectrums += [spectrum]
    return spectrums


def spectra_from_data(data):  # should be class function of Spectra
    spectrums = spectrums_from_data(data)
    print("spectra_from_daa: spectrums = " + str(spectrums))  # debugging
    return Spectra(spectrums=spectrums)


class Spectra:
    def __init__(
        self,
        file_path=None,
        folder=None,
        spectrums=None,
        data=None,
        tstamp=None,
        data_type="PVMS",
    ):
        #print(spectrums)  # debugging
        if file_path is not None and file_path[-4:] == ".pkl":
            with open(file_path, "rb") as f:
                data = pickle.load(f)
        elif data is None and spectrums is None:
            if data_type == "PVMS":
                data = read_PVMS(file_path)
                spectrums = data_to_spectrums(data)
            else:
                print(
                    "Spectra.__init__ does not yet support reading spectrums "
                    + "from files with data_type = "
                    + data_type
                )
        self.spectrums = spectrums
        self.data = data
        self.x = spectrums[0].x
        self.spectra = np.stack([spectrum.y for spectrum in spectrums])
        self.file_path = file_path
        if file_path is not None:
            self.folder, self.file = os.path.split(file_path)

    def __getitem__(self, key):
        if type(key) is int:
            return self.spectrums[key]
        elif key in self.data:
            return self.data[key]
        elif hasattr(self, key):
            return getattr(self, key)
        raise KeyError(
            "Spectra has no attribute " + key + ", and spectra.data has no key " + key
        )

    def save(self, file_name):
        data = {"x": self.x, "spectra": self.spectra}
        with open(file_name, "wb") as f:
            pickle.dump(data, f)


def data_to_spectrums(data):
    spectrums = []
    x_list = [col for col in data.keys() if re.search("^" + float_match + "$", col)]
    x = np.array([eval(v) for v in x_list])
    I_sort = np.argsort(x)
    x = x[I_sort]
    tstamp = data["tstamp"]
    ts = data[data["t_str"]]
    for (n, t) in enumerate(ts):
        try:
            y = np.array([float(data[x_list[I]][n]) for I in I_sort])
        except TypeError:
            print(
                "Warning!!! couldn't convert value "
                + " in spectrum number "
                + str(n)
                + ". Skipping that spectrum."
            )
            continue
        spectrum = Spectrum(x=x, y=y, t=t, tstamp=tstamp)
        spectrums += [spectrum]
    return spectrums
