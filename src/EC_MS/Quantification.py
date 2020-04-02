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


from .Molecules import Molecule
from .Combining import cut, get_cols_for_mass
from .Object_Files import lines_to_dictionary, date_scott
from .EC import sync_metadata
from . import Chem

molecule_directory = (
    os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "molecules"
)
data_directory = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data"
preferencedir = os.path.dirname(os.path.realpath(__file__)) + os.sep + "preferences"

with open(preferencedir + os.sep + "standard_colors.txt", "r") as f:
    lines = f.readlines()
    standard_colors = lines_to_dictionary(lines)["standard colors"]


def rewrite_spectra(
    NIST_file="default",
    RSF_file="default",
    mols="all",
    writesigma=False,
    overwrite=False,
    writespectra=False,
    writeRSF=False,
):
    """
    Reformats NIST data copied to a text file given by 'file', and writes it to
    the molecule files in data_directory
    """
    if writesigma or writespectra:
        if NIST_file == "default":
            NIST_file = data_directory + os.sep + "NIST_spectra_data.txt"
        with open(NIST_file, "r") as f:
            lines = f.readlines()
        structure = lines_to_dictionary(lines)
        sigma_dict = structure[
            "sigma_100eV"
        ]  # dictionary with electron impact ionizations in Ang^2
        spectra_dict = structure[
            "Spectra"
        ]  # dictionary with QMS spectra as mass1,rel_val1 mass2,relval2
    if writesigma:
        for (mol, sigma) in sigma_dict.items():
            if not (mols == "all" or mol in mols):
                continue
            try:
                m = Molecule(mol)
            except FileNotFoundError:
                with open(molecule_directory + os.sep + mol + ".txt", "w") as f:
                    f.write("name: " + mol)
            if hasattr(m, "sigma") and m.sigma is not None:
                if not overwrite:
                    continue
            l = ("sigma_100eV", sigma)
            m.write(l)
    if writespectra:
        for (mol, specline) in spectra_dict.items():
            if not (mols == "all" or mol in mols):
                continue
            m = Molecule(mol)
            if (
                hasattr(m, "spectrum")
                and m.spectrum is not None
                and len(m.spectrum) > 0
            ):
                if not overwrite:
                    continue
            masses = []
            rel_vals = []
            date = date_scott()
            spectrum = {"Title": "NIST_" + date}
            for spec in specline.split(" "):
                (mz, v) = spec.split(",")
                masses += ["M" + mz]
                rel_vals += [eval(v)]
            rel_vals = np.array(rel_vals)
            rel_vals = (
                100 * rel_vals / np.max(rel_vals)
            )  # normalize spectrum to % of max peak
            for (mass, rel_val) in zip(masses, rel_vals):
                spectrum[mass] = rel_val
            l = ("Spectrum", spectrum)
            m.write(l)
    if writeRSF:
        if RSF_file == "default":
            RSF_file = data_directory + os.sep + "Relative_Sensativity_Factors.txt"
        with open(RSF_file) as f:
            lines = f.readlines()
        structure = lines_to_dictionary(lines)
        for (mol, value) in structure.items():
            if not (mols == "all" or mol in mols):
                continue
            m = Molecule(mol)
            mass = "M" + str(value[0])
            RSF = value[1]
            l = ("Hiden", (mass, RSF))
            m.write(l)
    return structure


def calibration_compare(calmol=["H2", "O2", "CO2", "Cl2"]):
    """
    Checks the calibration factors calculated internally at the sniffer setup
    against predicted sensitivities based on the NIST database.
    There is no correlation. We are crazy sensitive to H2.
    """
    for mol in calmol:
        m = Molecule(mol, verbose=False)
        m.get_RSF()
        m.calibration_fit(verbose=False, ax=None)
        p = m.primary
        c1 = m.F_cal
        c2 = m.ifcs

        r12 = c1 / c2
        RSFit = "rsf" in dir(m)
        if RSFit:
            c3 = m.rsf
            r13 = c1 / c3

        print(
            mol
            + " at "
            + p
            + ", calibration factors:\n"
            + "\tc1 = "
            + str(c1)
            + " C/mol (experimental from Pt) \n"
            + "\tc2 = "
            + str(c2)
            + " Ang^2 (calculated from NIST) \n"
            + "r12 = "
            + str(r12)
            + "  (ratio) "
        )
        if RSFit:
            print(
                "\tc3 = "
                + str(c3)
                + "   (relative sensitivity factor from Hiden) \n"
                + "r13 = "
                + str(r13)
                + "  (ratio)"
            )


def RSF_to_F_cal(*args, **kwargs):
    print("'RSF_to_F_cal' has been renamed 'recalibrate'. Remember that next time!")
    return recalibrate(*args, **kwargs)


def recalibrate(
    quantify={},  # molecules we want to calc F_cal (at the given mass) for by extrapolation
    trust=None,  # says what to trust, e.g., 'internal', 'external', or 'all'
    trusted=[],  # list of Molecule objects with trusted F_cal
    internal=[],  # list of Molecule objects with F_cal from internal calibration
    external=[],  # list of Molecule objects with F_cal from internal calibration
    RSF_source="NIST",  #'NIST' uses partial ionization cross section
    # used in Trimarco2018 (Sniffer 2.0)
    #'Hiden' is RSF's from Hiden Analytical
    transmission_function="default",  # 'default' is T(M) = M^(-1/2)
    trendline=True,
    ax="new",
    labels=False,
    writeit=False,
    writeprimary=False,
    rewriteit=False,
):
    """
    ---- You need -----
    Calibration factor for a given set of molecules at a given set
    of masses.
    This function returns mdict, which is a dictionairy of the form
    {name:molecule, ...} where molecules are objects of the class EC_MS.Molecule,
    molecule.primary is the mass at which it is calibrated, and molecule.F_cal
    is its calibration factor at that mass in C/mol


    ------- You have -------
    - The names of the molecules you need and the masses you would like them
    Input these as a dictionairy:
        quantify = {molecule1:mass1, molecule2:mass2, ...}, e.g.:
        quantify = {'CH4':'M15', 'C2H4':'M26'}

    - Molecules for which you have a trusted calibration. Input these
    as a list of objects of the class EC_MS.Molecule. This list can be input
    as "trusted", "internal" (for internal calibrations), or "external". If
    both "internal" and "external" are given, the function by default only trusts
    the internal calibrations when calculating calibration factors for molecules
    in quantify, but co-plots both of them. e.g.:
        internal = {H2, CO2},
        external = {He, CO}
        where, for example, H2.primary = 'M2', H2.F_cal = 1.72, ...

    ---- additional options ----
    RSF_source: by default, uses ionization cross section and spectra from NIST,
        and a transmission function
    transmission_function: proportaional to probability an ion will make it
        through quadruopole. By default 1/sqrt(M)
    ax: by default, makes a new axis to plot calibration factors vs relative
        sensitivity factors
    trendline: whether to draw a trendline on the plot
    writeit: whether to save the calibration factors to the molecules' data files

    """
    print("\n\nfunction 'recalibrate' at your service!\n")

    # ----------- parse inputs, prepare output -------- #

    # prepare mdict, put trusted stuff in it
    mdict = {}  # to be returned by the function
    if quantify in ["all", "standard"]:
        quantify = [
            ("H2", "M2"),
            ("He", "M4"),
            ("CH4", "M15"),
            ("H2O", "M18"),
            ("N2", "M28"),
            ("CO", "M28"),
            ("C2H4", "M26"),
            ("C2H6", "M30"),
            ("O2", "M32"),
            ("Ar", "M40"),
            ("CO2", "M44"),
            ("Cl2", "M70"),
        ]

    if type(quantify) is list:
        quantify = dict(quantify)
    quantmols = list(quantify.keys())
    for name, mass in quantify.items():
        mdict[name] = Molecule(name)
        mdict[name].primary = mass
    if type(internal) is dict:
        print("converting internal from dict to list")
        internal = list(internal.values())
    if type(external) is dict:
        print("converting external from dict to list")
        external = list(external.values())
    if type(trusted) is dict:
        print("converting external from dict to list")
        trusted = list(trusted.values())
    if len(trusted) == 0 and trust is None:
        trusted = internal
    for m in external + internal + trusted:
        try:
            name = m.name
        except AttributeError:
            print("WARNING! " + str(m) + " has no name")
            pass
        else:
            if name in mdict:
                print(
                    "WARNING! recalibrate recieved multiple input molecules named "
                    + name
                )
            mdict[m.name] = m

    # store the mass and calibration factor of the trusted molecule in calmol and F_cals, respectively
    # (this could perhaps be done smarter)
    if trust == "all" or trust == "internal":
        trusted += internal
    if trust == "all" or trust == "external":
        trusted += external
    if trust == "files":  # trusts F_cal saved in molecule files. never used.
        trusted = [Molecule(name) for name in trusted]
    if len(internal) == 0 and len(external) == 0 and len(trusted) > 0:
        internal = trusted  # so that trusted points show up as squares
    calmol = {}
    F_cals = {}
    for m in trusted:
        try:
            calmol[m.name] = m.primary
            F_cals[m.name] = m.F_cal
        except AttributeError:
            print(
                "Cannot use "
                + str(m)
                + " to calibrate. Calibration must"
                + " be based on molecule objects with attributes 'primary' and 'F_cal'."
            )
    print("trusted calibrations at: " + str(calmol))
    print("F_cals = " + str(F_cals))
    if len(internal) == 0 and len(external) == 0:
        internal = calmol  # so that they plot as squares

    # --------- get the F_cal vs RSF relationship for the trusted molecules ------- #
    RSF_vec = []
    F_cal_vec = []

    if transmission_function == "default":

        def T(M):
            return M ** (-1 / 2)

    elif transmission_function == 1:

        def T(M):
            return 1

    else:
        T = transmission_function

    for m in trusted:
        name, mass, F_cal = m.name, m.primary, m.F_cal
        rsf = m.get_RSF(RSF_source=RSF_source, transmission_function=T, mass=mass)
        print(
            "F_"
            + name
            + "_"
            + mass
            + " = "
            + str(F_cal)
            + ", "
            + "rsf_"
            + name
            + "_"
            + mass
            + " = "
            + str(rsf)
        )
        F_cal_vec += [F_cal]
        RSF_vec += [rsf]
        if writeit:
            m.write(
                "#the following F_cal value is for "
                + mass
                + ", trusted on "
                + date_scott()
            )
            l = ("F_cal", F_cal)
            m.write(l)

    def fit_fun(x, a):
        return a * x

    if len(trusted) <= 1:
        r = F_cal_vec[0] / RSF_vec[0]
    else:
        r, pcov = curve_fit(fit_fun, RSF_vec, F_cal_vec, p0=1)
        try:
            r = r[0]  # I think it comes back as an array, but I just want a number
        except TypeError:
            pass
    RSF_unit = "a.u."  # {'Hiden':'a.u.', 'NIST':'a.u.'}[RSF_source]
    print(
        "\n--- Calibration Factor / rsf = " + str(r) + " (C/mol)/" + RSF_unit + "---\n"
    )

    # ----------- prepare the figure, plot the given F_cals
    if ax == "new":
        fig, ax = plt.subplots()
    if ax is not None:
        ax.set_xlabel("Relative Sensitivity Factor / [" + RSF_unit + "]")
        ax.set_ylabel("F_cal / [C/mol]")
        for m in internal:
            name, mass, F_cal = m.name, m.primary, m.F_cal
            rsf = m.get_RSF(
                RSF_source=RSF_source, transmission_function=T, mass=mass, verbose=False
            )
            try:
                color = m.get_color()
            except AttributeError:
                color = standard_colors[mass]
            print("plotting " + name + " as a color=" + color + " square.\n")
            ax.plot(rsf, F_cal, "s", color=color, markersize=10)
            if labels:
                ax.annotate(
                    name + " at m/z=" + mass[1:], xy=[rsf + 0.05, F_cal], color=color
                )
        for m in external:
            name, mass, F_cal = m.name, m.primary, m.F_cal
            rsf = m.get_RSF(
                RSF_source=RSF_source, transmission_function=T, mass=mass, verbose=False
            )
            try:
                color = m.get_color()
            except AttributeError:
                color = standard_colors[mass]
            print("plotting " + name + " as a color=" + color + " triangle.\n")
            ax.plot(rsf, F_cal, "^", color=color, markersize=10)
            if labels:
                ax.annotate(
                    name + " at m/z=" + mass[1:], xy=[rsf + 0.05, F_cal], color=color
                )

    # -------- use rsf to predict F_cal for all the other molecules
    rsf_max = 0
    for (name, m) in mdict.items():
        print("\n\n --- working on " + name + " ----")
        if name in quantify:
            mass = quantify[name]
        else:
            mass = m.primary
        color = m.get_color()

        # calculate RSF
        rsf = m.get_RSF(RSF_source=RSF_source, transmission_function=T, mass=mass)
        rsf_max = max(rsf_max, rsf)
        if writeit:
            m.write(
                "#the folowing rsf is calculated for " + mass + " on " + date_scott()
            )
            m.write(("rsf", rsf))
        if rsf is None:
            print("missing rsf for " + name + " at " + mass)
            continue  # then nothing to do

        # get (already plotted) or calculate (and plot) F_cal
        if name in F_cals:
            F_cal = F_cals[name]
        elif name in quantmols:
            F_cal = r * rsf  # the extrapolation!

        if name in quantmols and ax is not None:
            print("plotting " + name + " as a color=" + color + " dot.")
            ax.plot(rsf, r * rsf, ".", color=color, markersize=10)
            if labels:
                ax.annotate(
                    name + " at m/z=" + mass[1:], xy=[rsf + 0.05, F_cal], color=color
                )
            print(name + ": F_cal = " + str(F_cal))

            # write it to the Molecule object
            if name in quantmols:  # only do it for the molecules that were asked for
                if "cal" in m.__dict__:
                    m.cal[mass] = F_cal
                if "primary" not in m.__dict__:
                    m.primary = mass
                if m.primary == mass:
                    m.F_cal = F_cal

        # write it to the Molecule's data file
        if writeit:
            m.write(
                "#the following F_cal value is for "
                + mass
                + ", extrapolated "
                + "from trusted values based on RSF from "
                + RSF_source
                + " on "
                + date_scott()
            )
            l = ("F_cal", F_cal)
            m.write(l)
        if writeprimary:
            l = ("primary", mass)
            m.write(l)
        if rewriteit:
            m.rewrite()

    # make a trendline (needs to come here for use of rsf_max)
    if ax is not None and trendline:
        ax.plot([0, rsf_max], [0, r * rsf_max], "k--")

    # ---- done!
    print("\nfunction 'recalibrate' finished!\n\n")
    if ax is not None:
        return mdict, ax
    return mdict


def line_through_zero(x, y):
    """
    This returns a minimizing the square error of y = a * x
    """
    a_i = np.mean(y / x)
    pars_i = [a_i]

    def ax(x, a):
        return a * x

    pars, pcov = curve_fit(ax, x, y, p0=pars_i)
    #    pars = [tau, y0, y1]
    a = pars[0]
    return a


def get_signal(
    MS_data,
    mass,
    tspan=None,
    removebackground=None,
    background=None,
    t_bg=None,
    endpoints=5,
    fillcolor=None,
    unit="A",
    verbose=True,
    override=False,
    plotit=False,
    ax="new",
    return_bg=False,
):
    """
    Returns [x, y] where x is the time and y is QMS signal.
    A bit trivial, but I like having this function to work in parrallel
    to get_flux.
    """

    if verbose:
        print("geting signal for " + mass)

    xcol, ycol = get_cols_for_mass(mass, MS_data)
    x = MS_data[xcol]
    y = MS_data[ycol]

    if len(x) == 0:
        print("WARNIGN: no data for " + mass)
        return x, y

    if unit[-1] == "A":
        if unit[:-1] == "n" or unit[:-1] == "nano":
            y = y * 1e9
        elif unit[:-1] == "u" or unit[:-1] == "micro":
            y = y * 1e6
        elif unit[:-1] == "p" or unit[:-1] == "pico":
            y = y * 1e12

    if tspan is None:
        tspan = "tspan"
    if type(tspan) is str and not tspan == "all":
        try:
            tspan = MS_data[tspan]
        except KeyError:
            print("WARNING: no tspan available to get_signal()! using tspan='all'")
            tspan = "all"
    if not tspan == "all":
        try:
            x, y = cut(x, y, tspan, override=override)
        except TypeError:
            # print('x = ' + str(x) + ', y = ' + str(y) + ', tspan = ' + str(tspan)) # debugging
            x = tspan
            try:
                y = np.interp(tspan, x, y)
            except ValueError:
                print("WARNING: couldn't cut according to tspan=" + str(tspan))

    if len(x) == 0:
        print("WARNING: no signal in the requested tspan for " + mass)
        return [x, y]

    if background is None and t_bg is not None:
        background = "constant"
    if removebackground is None:
        removebackground = not (background is None)

    # print('background = ' + str(background)) # debugging
    if removebackground:
        if background is None:
            background = "constant"
        if background == "start":
            background = np.average(y[:endpoints])
        elif background == "constant":
            if type(removebackground) is float:
                background = removebackground * min(y)
            elif t_bg is not None:
                try:
                    if verbose:
                        print("Averaging background at t_bg = " + str(t_bg))
                    # mask = np.logical_and(t_bg[0]<x, x<t_bg[-1])
                    x_bg, y_bg = get_signal(
                        MS_data,
                        mass=mass,
                        tspan=t_bg,
                        removebackground=False,
                        unit=unit,
                    )
                    background = np.mean(y_bg)
                except TypeError:
                    if verbose:
                        print("Interpolating background at t_bg = " + str(t_bg))
                    background = np.interp(t_bg, x, y)
            else:
                background = min(y)
        elif type(background) is float:
            # background = background
            pass
        elif background == "linear":
            x_end = [np.average(x[:endpoints]), np.average(x[-endpoints:])]
            y_end = [np.average(y[:endpoints]), np.average(y[-endpoints:])]
            background = np.interp(x, x_end, y_end)

    # print('background = ' + str(background)) # debugging
    if plotit:
        if ax == "new":
            fig, ax = plt.subplots()
        ax.plot(x, y, "k")
        if removebackground:
            ax.plot(x, background * np.ones(x.shape), "r--")
            if fillcolor:
                ax.fill_between(x, background, y, where=y > background, color=fillcolor)
        # ax.set_title(mass)

    if removebackground:
        # y = y - 0.99 * background #so that we don't break the log scale.
        y = y - background
        # I should get rid of this and assume the caller knows what they're doing.

    if return_bg:
        return [x, y, background]
    else:
        return [x, y]


def get_flux(MS_data, mol, **kwargs):
    """
    returns [x, y] where x is the t corresponding to the primary mass of the
    molecule in 'mol' and y is the molecular flux in nmol/s, calculated from
    the MS_data for its primary mass and the value of F_cal read from the
    molecule's text file.

    Now moved to inside the class Molecules.Molecule
    """

    if type(mol) is str:
        m = Molecule(mol, verbose=False)
    else:
        m = mol

    return m.get_flux(MS_data, **kwargs)


def get_current(EC_data, tspan="tspan", unit="A", verbose=False):
    """
    Returns current in requested unit (default is A) over a requested tspan.
    I'm not happy with this function. I need to completely upgrade the way
    mass-, area-, and otherwise normalized currents and units are handled.
    This should work for now.
    """
    if verbose:
        print("getting current in " + unit)
    t_str = "time/s"
    V_str, J_str = sync_metadata(EC_data, verbose=False)
    if "/" in unit:  # then it's normalized in some way
        t, j = EC_data[t_str], EC_data[J_str]
    else:  # Then you want the raw current data
        try:
            I_str = EC_data["I_str"]
        except KeyError:
            I_str = "I/mA"
        t, j = EC_data[t_str], EC_data[I_str]
        if unit == "A":
            j = j * 1e-3

    if unit == "A/m^2":  # SI units
        j = j * 10  # mA/cm^2 --> A/m^2

    print(
        "Scott has not yet done something smart for current units. If you try"
        + " anything unusual it will likely mess up!"
    )

    if type(tspan) is str and not tspan == "all":
        tspan = EC_data[tspan]
    if not tspan == "all":
        t, j = cut(t, j, tspan)

    return [t, j]


def get_potential(EC_data, tspan="tspan", scale="RHE", verbose=False):
    """
    I'm not happy with this function. I need to completely upgrade the way
    mass-, area-, and otherwise normalized currents and units are handled.
    This should work for now.
    """
    if verbose:
        print("getting potential on " + scale + " scale if possible")
    t_str = "time/s"
    V_str, J_str = sync_metadata(EC_data, verbose=False)
    if scale == "RHE":  # then it's normalized in some way
        t, V = EC_data[t_str], EC_data[V_str]
    elif scale == "RE":  # Then you want the raw current data
        t, V = EC_data[t_str], EC_data[EC_data["E_str"]]

    print(
        "Scott has not yet done something smart for current units. If you try"
        + " anything unusual it will likely mess up!"
    )

    if type(tspan) is str and not tspan == "all":
        tspan = EC_data[tspan]
    if not tspan == "all":
        t, V = cut(t, V, tspan)

    return [t, V]


def predict_current(
    EC_and_MS,
    mols,
    tspan=None,
    RE_vs_RHE=None,
    A_el=None,
    ax="new",
    colors=None,
    verbose=1,
):
    """
    calculates a predicted electrical current based on MS_data and molecular
    data loaded from the files named in the list 'mols.'
    As of 16K28 just uses the primary masses
    tspan = x1 uses the time axis of the first primary mass.
    """
    V_str, J_str = sync_metadata(EC_and_MS, RE_vs_RHE, A_el)
    A_el = EC_and_MS["A_el"]
    if A_el is None:
        A_el = 1

    partials = []
    for mol in mols:
        molecule = Molecule(mol)
        """
        mass = molecule.primary
        x = EC_and_MS[mass + '-x']
        y = EC_and_MS[mass + '-y']
        F_cal = molecule.F_cal
        i_mol = y / F_cal       #molecular flux in mol/s
        """
        [x, i_mol] = get_flux(EC_and_MS, molecule)
        j_mol = i_mol / A_el  # molecular flux densith in mol/cm^2/s
        j = j_mol * molecule.n_el * Chem.Far * 1e3  # current density in mA/cm^2
        partials += [[x, j]]

    # this handling of tspan might be excessive for now...
    # it's intricate because I might need to interpolate MS data to EC time
    if tspan == None and "tspan_2" not in EC_and_MS.keys():
        tspan == "x1:"
    if tspan == "x1":
        t = partials[0][0]
    elif "time/s" in EC_and_MS.keys():
        t = EC_and_MS["time/s"]
    if tspan == None or tspan == "tspan_2":
        tspan = EC_and_MS["tspan_2"]
    if len(tspan) == 2 and type(tspan) is not str:
        t = [t_i for t_i in t if tspan[0] < t_i and t_i < tspan[1]]
    else:
        tspan = [t[0], t[-1]]

    if ax == "new":
        fig = plt.figure()
        ax = fig.add_subplot(111)

    js = []
    j_total = np.zeros(np.shape(t))
    for ([x, j], mol) in zip(partials, mols):
        I_keep = [I for (I, x_I) in enumerate(x) if t[0] < x_I and x_I < t[-1]]
        x = x[I_keep]
        j = j[I_keep]
        j_int = np.interp(t, x, j)
        j_total += j_int
        js += [j_int]
        if ax is not None and colors is not None:
            ax.plot(t, j_int, colors[mol], label=mol)

    if ax is not None:
        ax.plot(t, j_total, "k:", label="total")


if __name__ == "__main__":
    plt.close("all")
    # calibration_compare()

    from EC_MS import set_figparams

    set_figparams(figwidth=8)

    F_cals = {
        "CO2": 23.371483271036237,
        "H2": 29.899494846007542,
        "O2": 14.756572997297784,
    }
    mdict, ax = RSF_to_F_cal(
        calmol={"CO2": "M44"},
        RSF_source="Hiden",
        trust=["H2", "O2", "CO2"],
        F_cals=F_cals,
    )
    # calmol is the molecule(s) who's calibration is extrapolated to the others
    # RSF_source can be 'NIST' or 'Hiden'.
    # trust gives the molecules that become the red dots.
    # F_cals overrides the stored calibration factors for given molecules.
    # mdict containes the Molecule objects with the extrapolated F_cal
    # ax is the axis it's plotted on.

    for mol, m in mdict.items():
        rsf = m.get_RSF(RSF_source="Hiden", verbose=False)
        print(mol + "\tF_cal = " + str(m.F_cal) + "\tRSF = " + str(rsf))

    ax.set_xlabel("Relative sensitivity factor / [a.u.]")
    ax.set_ylabel("F$_{\mathrm{cal}}$ / [C mol$^{-1}$]")

    plt.savefig("Fcal_vs_RSF.png")
