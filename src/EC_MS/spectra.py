# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:55:31 2020

@author: scott
"""
import os, re, pickle
import numpy as np
from matplotlib import pyplot as plt

from .PVMassSpec import read_PVMS
from .MS import Peak
from .dataset import Dataset

float_match = "[-]?\d+[\.]?\d*(e[-+]?\d+)?"  # matches floats like '-3.5e4' or '7' or '245.13' or '1e-15'


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
        # print("Initiating Spectrum!")  # debugging
        if x is None and y is None:
            spectra = Spectra(file_path, data_type=data_type)
            spectrum = spectra[index]
            x, y, t, tstamp = spectrum.x, spectrum.y, spectrum.t, spectrum.tstamp
        # print(file_path)  # debugging
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

    def get_peak(self, Mspan=None, mass=None, width=None):
        if Mspan is None:
            M = float(mass[1:])
            Mspan = [M - width / 2, M + width / 2]
        x, y = self.get_signal(Mspan=Mspan)
        peak = Peak(x, y)
        return peak

    def __sub__(self, spectrum_2):
        y_diff = self.y - spectrum_2.y
        x = self.x
        t = self.t
        tstamp = self.tstamp
        return Spectrum(x=x, y=y_diff, t=t, tstamp=tstamp)


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
        name="spectra",
    ):
        # print(spectrums)  # debugging
        if file_path is not None and file_path[-4:] == ".pkl":
            with open(file_path, "rb") as f:
                spectra_data = pickle.load(f)
            tstamp = spectra_data["tstamp"]
            data = spectra_data["data"]
            self.x = spectra_data["x"]
            self.spectra = spectra_data["spectra"]
            if "name" in spectra_data:
                self.name = spectra_data["name"]
            self.data = data
            if tstamp is None and "tstamp" in data:
                tstamp = data["tstamp"]
            try:  # okay, doing this twice, but whatever.
                self.t = data[data["t_str"]]
            except KeyError:
                print("Warning!!! can't find t in self.data")
            spectrums = self.spectrums_from_spectra()
        elif data is None and spectrums is None:
            if data_type == "PVMS":
                data = read_PVMS(file_path)
                spectrums = data_to_spectrums(data)
                if tstamp is None and "tstamp" in data:
                    tstamp = data["tstamp"]
            else:
                print(
                    "Spectra.__init__ does not yet support reading spectrums "
                    + "from files with data_type = "
                    + data_type
                )
        self.tstamp = tstamp
        self.file_path = file_path
        if file_path is not None:
            self.folder, self.file = os.path.split(file_path)
            name = "spectra from " + self.file
        self.spectrums = spectrums
        self.data = data
        try:
            self.t = data[data["t_str"]]
        except KeyError:
            print("Warning!!! can't find t in self.data")
        if not hasattr(self, "x"):  # it'll already have this if loaded form pickle
            self.x = spectrums[0].x
        if not hasattr(
            self, "spectra"
        ):  # it'll already have this if loaded from pickle
            self.spectra = np.stack([spectrum.y for spectrum in spectrums])
        if not hasattr(self, name):
            self.name = name

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

    def __len__(self):
        return len(self.spectrums)

    def spectrums_from_spectra(self, spectra=None, x=None):
        if spectra is None:
            spectra = self.spectra
        if x is None:
            x = self.x
        spectrums = []
        for i, y in enumerate(spectra):
            if "t" in self.data:
                t_i = self.data["t"][i]
            elif hasattr(self, "t"):
                t_i = self.t[i]
            else:
                t_i = None
            # print(f"{t_i} should be {self.t[i]}") # debugging
            spectrum = Spectrum(x=x, y=y, t=t_i)
            spectrums += [spectrum]
        self.spectrums = spectrums
        return spectrums

    def save(self, file_name):
        """
        Has the problem that it won't save properly if spectrums are different
        length. However, PVMassSpec can't save that kind of file, and neither
        can (as of now) Zilien.
        """
        spectra_data = {
            "x": self.x,
            "spectra": self.spectra,
            "data": self.data,
            "tstamp": self.tstamp,
            "file_path": self.file_path,
            "name": self.name,
        }
        with open(file_name, "wb") as f:
            pickle.dump(spectra_data, f)

    def heat_plot(
        self,
        ax="new",
        vs="number",
        logscale=True,
        orientation="yx",
        zrange=None,
        **kwargs,
    ):
        """
        kwargs are passed on to imshow.
        """

        spectra = self.spectra
        if ax == "new":
            fig, ax = plt.subplots()

        if vs == "number":
            t = np.arange(len(self))
            t_label = "scan number"
        elif vs == "t":
            t = self.t
            t_label = "time / [s]"
        M = self.x

        # this makes the extent one increment off.
        t = np.append(t, 2 * t[-1] - t[-2])
        M = np.append(M, 2 * M[-1] - M[-2])

        if orientation == "xy":
            spectra = np.swapaxes(spectra, 0, 1)
            spectra = np.flip(spectra, axis=0)
            ax.set_ylabel("m/z")
            ax.set_xlabel(t_label)
        else:
            ax.set_ylabel(t_label)
            ax.set_xlabel("m/z")
        if "extent" not in kwargs:
            if orientation == "xy":
                extent = [t[0], t[-1], M[0], M[-1]]
            elif orientation == "yx":
                extent = [M[0], M[-1], t[-1], t[0]]
            kwargs.update(extent=extent)

        if logscale:
            spectra = np.log(spectra)
        if zrange is None:
            good = np.logical_and(~np.isnan(spectra), ~np.isinf(spectra))
            # print('spectra = \n' + str(spectra)) # debugging
            low = np.min(spectra[good])
            high = np.max(spectra[good])
        else:
            low = zrange[0]
            high = zrange[1]
        spectra[spectra < low] = low
        spectra[spectra > high] = high
        spectra[np.isnan(spectra)] = low
        spectra[np.isinf(spectra)] = low

        if "aspect" not in kwargs:
            kwargs.update(aspect="auto")
        elif kwargs["aspect"] == "square":
            if orientation == "xy":
                aspect = (t[-1] - t[0]) / (M[-1] - M[0])
            elif orientation == "yx":
                aspect = (M[-1] - M[0]) / (t[-1] - t[0])
            kwargs.update(aspect=aspect)
        if "cmap" not in kwargs:
            kwargs.update(cmap="inferno")

        ax.imshow(spectra, **kwargs)

    def get_dataset(
        self,
        masses,
        mode="gauss_height",
        fit_width=1,
        y_bg=None,
        bg_mode=None,
        endpoints=2,
    ):
        dataset = Dataset(data_type=mode)
        dataset.data = {"timecols": {}}
        for mass in masses:
            xcol, ycol = mass + "-x", mass + "-y"
            x = np.array([])
            y = np.array([])
            M = float(mass[1:])
            Mspan = [M - fit_width / 2, M + fit_width / 2]
            for spectrum in self.spectrums:
                peak = spectrum.get_peak(Mspan=Mspan)
                try:
                    peak.fit_gauss(y_bg=y_bg, bg_mode=bg_mode, endpoints=endpoints)
                    height = peak.height

                except (IndexError, AttributeError):
                    print(
                        f"Warning!!! can't fit data in the range for {mass} at t~{spectrum.t}. putting nan."
                    )
                    y = np.nan
                x = np.append(x, spectrum.t)
                y = np.append(y, height)
            dataset.add_data_col(xcol, x, col_type=mode)
            dataset.add_data_col(ycol, y, col_type=mode)
            dataset.timecols[ycol] = xcol
            if not "spectrum number" in dataset.data:
                # ^ this should get called for the first mass column
                n_vec = np.arange(len(x))
                dataset.add_data_col(col="spectrum number", value=n_vec, timecol="M4-x")
        dataset.data["tstamp"] = self.data["tstamp"]
        dataset.data["data_type"] = "spectra"
        dataset["title"] = self.name
        dataset.empty = False

        return dataset
