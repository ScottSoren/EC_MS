# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:22:35 2020

@author: scott
"""
import os, json
import numpy as np
from matplotlib import pyplot as plt

data_dir = "../pickles/"
extraction_dir = "../extractions/"

from siCalibration import Calibration
from EC_MS import Dataset

# show up in git!

class Extraction(Dataset):
    """
    Written 20E03 as a first set of ideas for how to implement the Experiment
    class that I proposed some days ago to Reshma.
    """

    def __init__(
        self,
        name=None,
        dataset=None,
        data_file=None,
        data_files=None,
        tspan_experiment=None,
        tspan_exchange=None,
        tspan_extraction=None,
        t_bg=None,
        calibration=None,
        calibration_file="20A25_sniffer.json",
        electrolyte="16O",
        film="18O",
        element="Pt",
        tspan_ratio=None,
        alpha=None,
    ):

        if name is None:
            name = f"extraction {element}{film} in {electrolyte}"
        if calibration is None and calibration_file is not None:
            calibration = Calibration.load(calibration_file)
        self.calibration = calibration
        self.mdict = get_EC_MS_mdict(calibration)
        # Later we will just directly use the siQuant mdict!
        if dataset is None:
            if data_file is not None:  # <-- load one data file
                if os.sep not in data_file:
                    data_file = os.path.join(data_dir, data_file)
                dataset = Dataset(data_file)
            elif data_files is not None:  # <-- synchronize multiple datasets!
                dataset = Dataset()
                for data_file in data_files:
                    if os.sep not in data_file:
                        data_file = os.path.join(data_dir, data_file)
                    dataset = dataset + Dataset(data_file)
        self.dataset = dataset
        self.data = dataset.data
        self.tspan_experiment = tspan_experiment
        self.tspan_exchange = tspan_exchange
        self.tspan_extraction = tspan_extraction
        self.tspan_ratio = tspan_ratio
        self.electrolyte = parse_isotope(electrolyte)
        self.film = parse_isotope(film)
        # ^ parse_isotope makes sure it's e.g. 18O instead of O18
        self.element = element
        self.t_bg = t_bg
        if t_bg is not None:
            self.set_background(t_bg=t_bg)
        self.alpha = alpha
        if alpha is None and tspan_ratio is not None:
            self.get_alpha(tspan=tspan_ratio, ax=None)
        self.n_ex = {}  # will store extraction results

    def as_dict(self):
        self_as_dict = {}
        # fmt: off
        self_as_dict.update(
            data_file=self.data_file, calibration_file=self.calibration_file,
            tspan_experiment=self.tspan_experiment, tspan_exchange=self.tspan_exchange,
            tspan_extraction=self.tspan_extraction, tspan_ratio=self.tspan_ratio,
            alpha=self.alpha, n_ex=self.n_ex, t_bg=self.t_bg,
            electrolyte=self.electrolyte, film=self.film, element=self.element
        )
        # fmt: on
        return self_as_dict

    def save(self, extraction_file=None):
        if extraction_file is None:
            extraction_file = self.name + ".json"
        if os.sep in extraction_file:
            path_to_file = extraction_file
        else:
            path_to_file = os.path.join(extraction_dir, extraction_file)
        self_as_dict = self.as_dict()
        with open(path_to_file, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def load(cls, extraction_file):
        path_to_file = os.path.join(extraction_dir, extraction_file)
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        return cls(**self_as_dict)

    def plot_experiment(self, *args, **kwargs):
        """
        Just adds to the plot_experiment of the Dataset class that it uses
        the calibration and tspan_experiment by default
        """
        to_plot = []
        to_plot_0 = None
        if len(args) > 0:
            to_plot_0 = args[0]
        elif "mols" in kwargs:
            to_plot_0 = kwargs["mols"]
        if to_plot_0 is not None:
            for i, thing in enumerate(to_plot_0):
                if type(thing) is str and thing in self.mdict:
                    to_plot += [self.mdict[thing]]
                else:
                    to_plot += [thing]
        elif "masses" not in kwargs:
            to_plot = self.mdict

        # print(f"to_plot_0 = {to_plot_0}") # debugging
        # print(f"to_plot = {to_plot}") # debugging

        if len(args) > 0:
            args[0] = to_plot
        elif "mols" in kwargs or "masses" not in kwargs:
            kwargs["mols"] = to_plot

        if len(args) < 2 and "tspan" not in kwargs:
            kwargs.update(tspan=self.tspan_experiment)

        return super(Extraction, self).plot_experiment(*args, **kwargs)

    def get_alpha(self, tspan=None, t_bg=None, simple=True, ax="new"):
        if tspan is None:
            tspan = self.tspan_ratio

        if ax == "new":
            tspan_plot = [2 * tspan[0] - tspan[-1], 2 * tspan[-1] - tspan[0]]
            ax = self.plot_signal(
                ["M32", "M34", "M36"], tspan=tspan_plot, unit="pA", t_bg=t_bg
            )
        linewidth = 4

        x32, y32 = self.get_signal("M32", tspan=tspan, t_bg=t_bg, unit="pA")
        x34, y34 = self.get_signal("M34", tspan=tspan, t_bg=t_bg, unit="pA")
        x36, y36 = self.get_signal("M36", tspan=tspan, t_bg=t_bg, unit="pA")

        Y32 = np.mean(y32)
        Y34 = np.mean(y34)
        Y36 = np.mean(y36)

        print(f"self.electrolyte = {self.electrolyte}")  # debugging
        if simple:
            if self.electrolyte == "18O":
                gamma = Y34 / Y36
                self.gamma = gamma
                alpha = gamma / (2 + gamma)
                if ax is not None:
                    ax.plot(x34, y34, "r", linewidth=linewidth)
                    ax.plot(x36, y36, "g", linewidth=linewidth)
            elif self.electrolyte == "16O":
                beta = Y34 / Y32
                self.beta = beta
                alpha = 2 / (2 + beta)
                if ax is not None:
                    ax.plot(x32, y32, "k", linewidth=linewidth)
                    ax.plot(x34, y34, "r", linewidth=linewidth)
            else:
                print(
                    f"Warining!!! Can't get ratio with self.electrolyte={self.electrolyte}"
                )
                alpha = None
        else:
            from scipy.optimize import minimize

            y_hat = np.array([Y32, Y34, Y36])
            y_hat = y_hat / np.sum(y_hat)

            # g is the H2(16)O / H2(18)O ratio, called gamma elsewhere
            def sqerror(alpha):
                return (
                    (alpha ** 2 - y_hat[0]) ** 2
                    + (2 * alpha * (1 - alpha) - y_hat[1]) ** 2
                    + ((1 - alpha) ** 2 - y_hat[2]) ** 2
                )

            res = minimize(sqerror, 0.5)
            alpha = res.x[0]
            if ax is not None:
                ax.plot(x32, y32, "k", linewidth=linewidth)
                ax.plot(x34, y34, "r", linewidth=linewidth)
                ax.plot(x36, y36, "g", linewidth=linewidth)
        self.alpha = alpha
        return alpha

    def get_majors_and_minors(
        self, mol="O2",
    ):
        if self.electrolyte == "18O":
            if mol == "O2":
                majors = [self.mdict["O2_M36"]]
                minors = [self.mdict["O2_M34"], self.mdict["O2_M32"]]
            elif mol == "CO2":
                majors = [self.mdict["CO2_M46"], self.mdict["CO2_M48"]]
                minors = [self.mdict["CO2_M44"]]
        elif self.electrolyte == "16O":
            if mol == "O2":
                majors = [self.mdict["O2_M32"]]
                minors = [self.mdict["O2_M34"], self.mdict["O2_M36"]]
            elif mol == "CO2":
                majors = [self.mdict["CO2_M44"]]
                minors = [self.mdict["CO2_M46"], self.mdict["CO2_M48"]]
        return majors, minors

    def get_ratio(self):
        alpha = self.alpha
        if self.electrolyte == "18O":
            ratio = 2 * alpha / (1 - alpha)
        elif self.electrolyte == "16O":
            ratio = 2 * (1 - alpha) / (alpha)
        return ratio

    def plot_exchange(self, mol="O2", tspan=None, t_bg=None, axes="new", unit="pmol/s"):
        if tspan is None or tspan == "experiment":
            tspan = self.tspan_experiment
        elif tspan == "exchange":
            tspan = self.tspan_exchange
        elif tspan == "extraction":
            tspan = self.tspan_extraction
        ratio = self.get_ratio()
        majors, minors = self.get_majors_and_minors(mol=mol)

        if axes == "new":
            axes = self.plot_experiment(
                mols=[minors, majors],
                tspan=tspan,
                t_bg=t_bg,
                logplot=False,
                ax=axes,
                unit=unit,
            )
        else:
            for molecule in minors:
                self.plot_flux(
                    molecule, ax=axes[0], tspan=tspan, unit=unit, logplot=False
                )
            for molecule in majors:
                self.plot_flux(
                    molecule, ax=axes[-1], tspan=tspan, unit=unit, logplot=False
                )

        if True:  # highlight the labeled lattice oxygen evolution
            x1, y1 = self.get_flux(majors[0], t_bg=t_bg, unit="pmol/s", tspan=tspan)
            x2, y2 = self.get_flux(minors[0], t_bg=t_bg, unit="pmol/s", tspan=tspan)
            y2_interp = np.interp(x1, x2, y2)
            color_1 = minors[0].get_color()
            color_2 = majors[0].get_color()
            axes[0].fill_between(
                x1,
                y1 * ratio,
                y2_interp,
                color=color_1,
                where=y2_interp > y1 * ratio,
                alpha=0.5,
            )
            axes[0].fill_between(
                x1,
                y1 * ratio,
                y2_interp,
                color=color_2,
                where=y2_interp < y1 * ratio,
                alpha=0.2,
                hatch="//",
            )
        axes[-1].set_ylim([l / ratio for l in axes[0].get_ylim()])

        return axes

    def plot_extraction_vs_potential(
        self, mol="O2", tspan=None, unit="pmol/s", ax="new", reverse=True
    ):
        if tspan == None:
            tspan = self.tspan_extraction
        ratio = self.get_ratio()
        majors, minors = self.get_majors_and_minors(mol=mol)
        print(f"majors={majors}, minors={minors}")  # debugging
        ax = self.plot_vs_potential(
            mols=[minors, majors], tspan=tspan, unit=unit, ax=ax, logplot=False
        )
        ax[-1].set_ylim([l / ratio for l in ax[0].get_ylim()])
        if reverse:
            for ax_i in ax:
                ax_i.invert_xaxis()
        return ax

    def quantify_extraction(self, mol="O2", tspan=None, unit="pmol/s"):
        if tspan == None:
            tspan = self.tspan_extraction
        ratio = self.get_ratio()
        majors, minors = self.get_majors_and_minors(mol=mol)
        Y = 0
        for molecule in minors:
            x, y = self.get_flux(molecule, tspan=tspan, unit=unit)
            Y += np.trapz(y, x)
        for molecule in majors:
            x, y = self.get_flux(molecule, tspan=tspan, unit=unit)
            Y -= np.trapz(y, x) * ratio
        self.n_ex[mol] = Y
        return Y


def get_EC_MS_mdict(calibration):
    from EC_MS import Molecule

    mdict = {}
    for mol, molecule in calibration.mdict.items():
        if mol in calibration.real_names:
            mol_name = calibration.real_names[mol]
        else:
            mol_name = mol
        mdict[mol] = Molecule(mol_name)
        mdict[mol].name = mol
        F = molecule.F
        mdict[mol].F = F
        key_max = None
        F_max = 0
        for key, val in F.items():
            if val > F_max:
                F_max = val
                key_max = key
        mdict[mol].F_cal = F_max
        mdict[mol].primary = key_max
    return mdict


def parse_isotope(isotope):
    if isotope in ["18O", "O18"]:
        return "18O"
    elif isotope in ["16O", "O16"]:
        return "16O"
    else:
        return isotope
