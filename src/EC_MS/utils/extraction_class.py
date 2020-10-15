# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:22:35 2020

@author: scott
"""
import os, json
import numpy as np

# from matplotlib import pyplot as plt

STANDARD_DATA_DIR = "../pickles/"
STANDARD_EXTRACTION_DIR = "../extractions/"

from EC_MS import Dataset


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
        data_dir=STANDARD_DATA_DIR,
        tspan_experiment=None,
        tspan_exchange=None,
        tspan_extraction=None,
        RE_vs_RHE=None,
        A_el=None,
        t_bg=None,
        calibration=None,
        mdict=None,
        calibration_file=None,  # "20A25_sniffer.json",
        electrolyte="16O",
        film="18O",
        element="Pt",
        tspan_ratio=None,
        alpha=None,
        n_ex=None,
    ):

        if name is None:
            name = f"extraction {element}{film} in {electrolyte}"
        self.name = name
        self.calibration_file = calibration_file
        if (not calibration) and calibration_file and not mdict:
            print(f"requested calibration_file = {calibration_file}")
            raise NotImplementedError(
                "proper calibration is not implemented in EC_MS. "
                "You have to import the calibration externally. "
                "An EC_MS mdict should also work.\n"
                "...The functionality should come soon to ixdat."
            )
        self.calibration = calibration
        if calibration and not mdict:
            mdict = get_EC_MS_mdict(calibration)
        self.mdict = mdict
        # Later we will just directly use the siQuant mdict!
        self.data_file = data_file
        self.data_files = data_files
        self.data_dir = data_dir
        if dataset is None:
            if data_file is not None:  # <-- load one data file
                if os.sep not in str(data_file):
                    path_to_data = os.path.join(data_dir, data_file)
                else:
                    path_to_data = data_file
                dataset = Dataset(path_to_data)
            elif data_files is not None:  # <-- synchronize multiple datasets!
                dataset = Dataset()
                for data_file in data_files:
                    if os.sep not in data_file:
                        path_to_data = os.path.join(data_dir, data_file)
                    dataset = dataset + Dataset(path_to_data)
        self.dataset = dataset
        self.data = dataset.data
        self.RE_vs_RHE = RE_vs_RHE
        self.A_el = A_el
        if RE_vs_RHE is not None or A_el is not None:
            self.normalize(RE_vs_RHE=RE_vs_RHE, A_el=A_el)
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
        if n_ex is None:
            n_ex = {}
        self.n_ex = n_ex  # will store extraction results

    def as_dict(self):
        self_as_dict = {}
        # fmt: off
        self_as_dict.update(
            name=self.name,
            data_dir=self.data_dir, data_file=self.data_file,
            data_files=self.data_files, calibration_file=self.calibration_file,
            RE_vs_RHE=self.RE_vs_RHE, A_el=self.A_el,
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
            path_to_file = os.path.join(STANDARD_EXTRACTION_DIR, extraction_file)
        self_as_dict = self.as_dict()
        with open(path_to_file, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def load(cls, extraction_file, **kwargs):
        path_to_file = os.path.join(STANDARD_EXTRACTION_DIR, extraction_file)
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    def plot_experiment(self, *args, **kwargs):
        """
        Just adds to the plot_experiment of the Dataset class that it uses
        the calibration and tspan_experiment by default.
        Molecules can be given as str and the Extraction tries to look them
        up in its calibration.
        """
        to_plot = []  # this will be either mols or masses to plot
        to_plot_0 = None  # first we see if the function caller told us what to plot
        if len(args) > 0:
            to_plot_0 = args[0]  # first positional argument should be mols or masses!
        elif "mols" in kwargs:  # ... but mols can also be given as a kwarg
            to_plot_0 = kwargs["mols"]
        elif "masses" in kwargs:  # ... and so can masses
            pass  # but then we don't actually need to do anything to it

        if to_plot_0:
            # now we go through and put the requested calibrated objects in to_plot
            for i, thing in enumerate(to_plot_0):
                print(f"thing to plot = {thing}")  # debugging
                if isinstance(thing, str) and thing in self.mdict:
                    # excellent! then we've got the name of a calibrated molecule.
                    to_plot += [self.mdict[thing]]
                elif isinstance(thing, list) or isinstance(thing, tuple):
                    # ^ this will be the case if they ask for mols left and mols right
                    to_plot += [[]]  # we need to mirror the list in to_plot
                    for subthing in thing:
                        print(f"subthing to plot = {subthing}")  # debugging
                        if isinstance(subthing, str) and subthing in self.mdict:
                            # excellent! then we've got the name of a calibrated molecule.
                            to_plot[-1] += [self.mdict[subthing]]
                        else:
                            # then we assume they know what they're doing.
                            to_plot[-1] += [subthing]
                else:
                    # then we assume they know what they're doing.
                    to_plot += [thing]
        elif "masses" not in kwargs:
            # okay, so this is actually the case if they don't ask for anything to plot
            #   and then by default we try to plot everything in the calibration
            to_plot = self.mdict

        # print(f"to_plot_0 = {to_plot_0}") # debugging
        # print(f"to_plot = {to_plot}") # debugging

        if len(args) > 0:
            # then we overwrite to_plot_0, which is args[0], with to_plot:
            args[0] = to_plot
        elif "mols" in kwargs or "masses" not in kwargs:
            # this is the case if they asked for mols or didn't ask for anything.
            #   and then we've generated the calibrated plotted list:
            kwargs["mols"] = to_plot

        if len(args) < 2 and "tspan" not in kwargs:
            # ah, yes, we use the Experiment tspan by default instead of the Dataset tspan:
            kwargs.update(tspan=self.tspan_experiment)

        # and now we're ready to call plot_experiment via Dataset.plot_experiment!

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
        """
        Get the majority and minority isotopes of mol produced in the electrolyte

        majors are ^{16}O2 and C^{16}O2 in 16O electrolyte,
        ^{18}O2 and (C^{16}O^{18}O and C^{18}O2) in 18O electrolyte,
        minors are the other isotopes.
        """
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
            ratio = 2 * (1 - alpha) / alpha
        return ratio

    def plot_exchange(
        self, mol="O2", tspan=None, t_bg=None, axes="new", unit=None, **kwargs
    ):
        if tspan is None or tspan == "experiment":
            tspan = self.tspan_experiment
        elif tspan == "exchange":
            tspan = self.tspan_exchange
        elif tspan == "extraction":
            tspan = self.tspan_extraction
        ratio = self.get_ratio()
        majors, minors = self.get_majors_and_minors(mol=mol)

        unit_right = kwargs.get("unit_right", unit)
        unit_left = kwargs.get("unit_left", unit)

        if axes == "new":
            axes = self.plot_experiment(
                mols=[minors, majors],
                tspan=tspan,
                t_bg=t_bg,
                logplot=False,
                ax=axes,
                unit=unit,
                **kwargs,
            )
        else:
            for molecule in minors:
                self.plot_flux(
                    molecule, ax=axes[0], tspan=tspan, unit=unit_left, logplot=False
                )
            for molecule in majors:
                self.plot_flux(
                    molecule, ax=axes[-1], tspan=tspan, unit=unit_right, logplot=False
                )

        unit_ratio = 1
        if (
            unit_left
            and unit_right
            and unit_left.startswith("p")
            and unit_right.startswith("n")
        ):
            unit_ratio *= 1e3
        axes_ratio = ratio * unit_ratio
        if True:  # highlight the labeled lattice oxygen evolution
            x1, y1 = self.get_flux(majors[0], t_bg=t_bg, unit=unit_right, tspan=tspan)
            x2, y2 = self.get_flux(minors[0], t_bg=t_bg, unit=unit_left, tspan=tspan)
            y2_interp = np.interp(x1, x2, y2)
            color_1 = minors[0].get_color()
            color_2 = majors[0].get_color()
            axes[0].fill_between(
                x1,
                y1 * axes_ratio,
                y2_interp,
                color=color_1,
                where=y2_interp > y1 * ratio,
                alpha=0.5,
            )
            axes[0].fill_between(
                x1,
                y1 * axes_ratio,
                y2_interp,
                color=color_2,
                where=y2_interp < y1 * ratio,
                alpha=0.2,
                hatch="//",
            )
        axes[-1].set_ylim([l / axes_ratio for l in axes[0].get_ylim()])

        return axes

    def plot_extraction_vs_potential(
        self, mol="CO2", tspan=None, unit=None, ax="new", reverse=True
    ):
        if tspan is None:
            tspan = self.tspan_extraction
        ratio = self.get_ratio()
        majors, minors = self.get_majors_and_minors(mol=mol)
        print(f"majors={majors}, minors={minors}")  # debugging
        ax = self.plot_vs_potential(
            mols=[minors, majors], tspan=tspan, unit=unit, ax=ax, logplot=False
        )
        ax[-1].set_ylim([l / ratio for l in ax[0].get_ylim()])
        if reverse:
            ax[0].invert_xaxis()
            ax[1].invert_xaxis()
        return ax

    def create_excess_mol(self, mol="CO2", ratio=None, ratio_type=None, name=None):
        """Return EC_MS.Molecule object who's get_flux calculates lattice oxygen ex."""
        if name is None:
            name = "excess_" + mol
        majors, minors = self.get_majors_and_minors(mol=mol)
        if ratio is None:
            if ratio_type is None:
                if mol == "CO2" and self.electrolyte == "18O":
                    ratio_type = "single"
                else:
                    ratio_type = "ratio"
            if ratio_type == "single":  # expect only one O atom from electrolyte
                ratio = self.alpha / (1 - self.alpha)
            elif ratio_type == "alpha":  # expect only one O atom from electrolyte
                ratio = self.alpha
            else:
                ratio = self.get_ratio()
        print(f"ratio = {ratio}")  # debugging
        excess_molecule = minors[0]
        excess_molecule.name = name
        excess_molecule.cal_mat = {
            minors[0].primary: 1 / minors[0].F_cal,
            majors[0].primary: -1 / majors[0].F_cal * ratio,
        }
        self.mdict[name] = excess_molecule
        return excess_molecule

    def get_excess(self, mol="O2", x=None, tspan=None, unit="pmol/s"):
        majors, minors = self.get_majors_and_minors(mol=mol)
        x1, y1 = self.get_flux(minors[0], tspan=tspan, unit=unit)
        x0, y0 = self.get_flux(majors[-1], tspan=tspan, unit=unit)
        ratio = self.get_ratio()
        if x is None:
            if tspan is None:
                tspan = self.tspan_exchange
            x = np.linspace(tspan[0], tspan[-1], 100)
        y = np.interp(x, x1, y1) - np.interp(x, x0, y0) * ratio
        return x, y

    def quantify_extraction(self, mol="O2", tspan=None, unit="pmol/s"):
        x, y = self.get_excess(mol=mol, tspan=tspan, unit=unit)
        Y = np.trapz(y, x)
        self.n_ex[mol] = Y
        return Y


def get_EC_MS_mdict(calibration, mols=None, get_cal_mat=True):
    from EC_MS import Molecule

    mdict = {}
    if mols is None:
        mols = set([])
        if hasattr(calibration, "mol_list"):
            mols = mols.union(calibration.mol_list)
        if hasattr(calibration, "mdict"):
            mols = mols.union(calibration.mdict)  # adds mdict's keys to the set mols
    for mol in mols:
        print(f"geting EC_MS molecule {mol} from calibration")
        molecule = calibration.molecule(mol, with_Q=get_cal_mat)
        if mol in calibration.real_names:
            mol_name = calibration.real_names[mol]
        else:
            mol_name = mol
        mdict[mol] = Molecule(mol_name)
        mdict[mol].name = mol
        F = molecule.F
        mdict[mol].F = F
        if get_cal_mat:
            try:
                mdict[mol].cal_mat = molecule.Q
            except AttributeError:
                print(f"Warning!!! couldn't get Q for molecule {mol}.")
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
    if isotope in ["18O", "O18", "18", 18]:
        return "18O"
    elif isotope in ["16O", "O16", "16", 16]:
        return "16O"
    else:
        return isotope
