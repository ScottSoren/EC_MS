# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 19:07:45 2016
Most recently edited: 17B21

@author: scott
"""

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
import os

# from mpl_toolkits.axes_grid1 import make_axes_locatable

from .EC import sync_metadata, select_cycles

from .Quantification import get_flux, get_signal
from .Object_Files import lines_to_dictionary
from .Molecules import Molecule

preferencedir = os.path.dirname(os.path.realpath(__file__)) + os.sep + "preferences"
with open(preferencedir + os.sep + "standard_colors.txt", "r") as f:
    lines = f.readlines()
    standard_colors = lines_to_dictionary(lines, removecomments=False)[
        "standard colors"
    ]


def get_standard_colors():
    return standard_colors


def colorax(ax, color, lr="right", xy="y"):
    ax.spines[lr].set_color(color)
    ax.tick_params(axis=xy, color=color)
    ax.tick_params(axis=xy, labelcolor=color)
    if xy == "y":
        ax.yaxis.label.set_color(color)
    if xy == "x":
        ax.xaxis.label.set_color(color)


def align_zero(ax, ax_ref, xy="y"):
    ylim0 = ax.get_ylim()
    ylim_ref = ax_ref.get_ylim()
    A = ylim_ref[-1] / ylim_ref[0]
    B = ylim0[-1] - ylim0[0]
    a = B / (A - 1)
    b = A * B / (A - 1)
    ylim = [a, b]
    ax.set_ylim(ylim)
    return ylim


def smooth(y, n_points):
    smoother = np.ones((n_points,)) / n_points
    y_smooth = np.convolve(y, smoother, mode="same")
    return y_smooth


def plot_EC_vs_t(
    data,
    t_str="time/s",
    J_str=None,
    V_str=None,
    V_color="k",
    J_color="r",
    ax="new",
    verbose=True,
    **kwargs
):
    if J_str is None or V_str is None:
        V_str_0, J_str_0 = sync_metadata(data)
        if J_str is None:
            J_str = J_str_0
        if V_str is None:
            V_str = V_str_0
    t, V, J = data[t_str], data[V_str], data[J_str]
    if ax == "new":
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    else:
        ax1 = ax[0]
        ax2 = ax[1]
    ax1.plot(t, V, color=V_color, **kwargs)
    ax2.plot(t, J, color=J_color, **kwargs)
    ax1.set_xlabel(t_str)
    ax1.set_ylabel(V_str)
    ax2.set_ylabel(J_str)
    return [ax1, ax2]


def plot_vs_potential(
    CV_and_MS_0,
    colors=None,
    tspan=None,
    RE_vs_RHE=None,
    A_el=None,
    cycles="all",
    ax1="new",
    ax2="new",
    ax=None,  # spec='k-',
    overlay=0,
    logplot=True,
    leg=False,
    verbose=True,
    removebackground=None,
    background=None,
    t_bg=None,
    endpoints=3,
    masses="all",
    masses_left=None,
    masses_right=None,
    mols=None,
    mols_left=None,
    mols_right=None,
    unit=None,
    smooth_points=0,
    emphasis="ms",
    hspace=0.1,
    left_space=0.15,
    right_space=0.95,  # for gridspec
    J_str=None,
    V_str=None,
    fig=None,
    spec={},
    t_str=None,
    **kwargs
):
    """
    This will plot current and select MS signals vs E_we, as is the
    convention for cyclic voltammagrams. added 16I29

    #there's a lot of code here that's identical to plot_experiment. Consider
    #having another function for e.g. processing these inputs.
    """
    if verbose:
        print("\n\nfunction 'plot_vs_potential' at your service!\n")
    if type(logplot) is not list:
        logplot = [logplot, False]
    if removebackground is None:
        removebackground = t_bg is not None

    spec.update(kwargs)  # extra arguments are passed to plt.plot

    # prepare axes. This is ridiculous, by the way.
    CV_and_MS = CV_and_MS_0.copy()  # 17C01
    if "data_type" in CV_and_MS and CV_and_MS["data_type"][0:2] == "EC":
        ax = plot_vs_potential_EC(
            data=CV_and_MS,
            tspan=tspan,
            RE_vs_RHE=RE_vs_RHE,
            A_el=A_el,
            cycles="all",
            ax=ax,  # spec='k-',
            J_str=J_str,
            V_str=V_str,
            t_str=t_str,
            fig=fig,
            spec=spec,
            verbose=verbose,
        )
        return ax

    if ax == "new":
        ax1 = "new"
        ax2 = "new"
    elif ax is not None:
        ax1 = ax[0]
        ax2 = ax[1]
    if ax1 != "new":
        figure1 = ax1.figure
    elif ax2 != "new":
        figure1 = ax2.figure
    else:
        if fig is None:
            figure1 = plt.figure()
        else:
            figure1 = fig
    if overlay:
        if ax1 == "new":
            ax1 = figure1.add_subplot(111)
        if ax2 == "new":
            ax2 = ax1.twinx()
    else:
        if ax1 == "new":
            if emphasis == "MS":
                gs = gridspec.GridSpec(3, 1)
                # gs.update(hspace=0.025)
                ax1 = plt.subplot(gs[0:2, 0])
                ax2 = plt.subplot(gs[2:3, 0])
            elif emphasis == "ms":
                gs = gridspec.GridSpec(5, 1)
                # gs.update(hspace=0.025)
                ax1 = plt.subplot(gs[0:3, 0])
                ax2 = plt.subplot(gs[3:5, 0])
            elif emphasis == "EC":
                gs = gridspec.GridSpec(3, 1)
                # gs.update(hspace=0.025)
                ax1 = plt.subplot(gs[0, 0])
                ax2 = plt.subplot(gs[1:3, 0])
            else:
                gs = gridspec.GridSpec(8, 1)
                # gs.update(hspace=0.025)
                ax1 = plt.subplot(gs[0:4, 0])
                ax2 = plt.subplot(gs[4:8, 0])
            gs.update(hspace=hspace, left=left_space, right=right_space)

    if type(logplot) is int:
        logplot = [logplot, logplot]
    if logplot[0]:
        ax1.set_yscale("log")
    if logplot[1]:
        ax2.set_yscale("log")

    # get EC data
    V_str, J_str = sync_metadata(
        CV_and_MS, RE_vs_RHE=RE_vs_RHE, A_el=A_el, V_str=V_str, J_str=J_str
    )

    V = CV_and_MS[V_str]
    J = CV_and_MS[J_str]
    if t_str is None:
        if "t_str" in CV_and_MS:
            t_str = CV_and_MS["t_str"]
        else:
            t_str = "time/s"

    # get time variable and plotting indexes
    t = CV_and_MS[t_str]
    if tspan is None:  # then use the whole range of overlap
        try:
            tspan = CV_and_MS["tspan"]
        except KeyError:
            tspan = [min(t), max(t)]

    mask = np.logical_and(tspan[0] < t, t < tspan[-1])

    if ax2 is not None:
        # plot EC-lab data
        ec_spec = spec.copy()
        if "color" not in ec_spec.keys():
            ec_spec["color"] = "k"
        # print('len(data[' + V_str + '] = ' + str(len(V))) # debugging
        # print('len(data[' + J_str + '] = ' + str(len(J))) # debugging
        # print('len(mask) = ' + str(len(mask))) # debugging

        t_plot, V_plot, J_plot = t[mask], V[mask], J[mask]
        ax2.plot(V_plot, J_plot, **ec_spec)
        # maybe I should use EC.plot_cycles to have different cycles be different colors. Or rewrite that code here.
        ax2.set_xlabel(V_str)
        ax2.set_ylabel(J_str)

    axes = [ax1, ax2]

    # ---- parse inputs for (evt. calibrated) mass spectrometer to plot ------
    try:
        if type(mols[0]) in (list, tuple):
            mols_left = mols[0]
            mols_right = mols[1]
            mols = None
    except (IndexError, TypeError):
        pass
    try:
        if type(masses[0]) in (list, tuple):
            masses_left = masses[0]
            masses_right = masses[1]
            masses = None
    except (IndexError, TypeError):
        pass
    #    print(masses)
    colors_right = None
    colors_left = None
    if ax1 is not None:  # option of skipping an axis added 17C01
        # check if we're going to plot signals or fluxes:
        quantified = False  # added 16L15
        if mols is not None:
            quantified = True
            colors_left = mols  # added 17H11
        elif mols_left or mols_right is not None:
            quantified = True
            colors_left = mols_left
            colors_right = mols_right
        elif masses is not None:
            #            print('masses specified')
            quantified = False
            colors_left = masses
        elif masses_left or masses_right is not None:
            quantified = False
            colors_left = masses_left
            colors_right = masses_right
        elif (
            (type(colors) is dict and list(colors.keys())[0][0] == "M")
            or (type(colors) is list and type(colors[0]) is str and colors[0][0] == "M")
            or (type(colors) is str and colors[0] == "M")
        ):
            if verbose:
                print("uncalibrated data to be plotted.")
            masses = colors
            colors_left = masses
        else:
            quantified = True
            mols = colors

        if not quantified and masses == "all":  # this is now the default!
            masses = [
                key[:-2]
                for key in CV_and_MS.keys()
                if key[0] == "M" and key[-2:] == "-y"
            ]
            colors_left = masses

        if (
            colors_left is not None
            and type(colors_left) is not list
            and type(colors_left) is not dict
        ):
            colors_left = [colors_left]
        if (
            colors_right is not None
            and type(colors_right) is not list
            and type(colors_right) is not dict
        ):
            colors_right = [colors_right]
        #        print(type(colors))

        if unit is None:
            if quantified:
                unit = "pmol/s"
            else:
                unit = "pA"

        # then do it.
        if colors_right is not None:
            if ax is None or len(ax) < 3:  # so I can reuse right axis
                axes += [axes[0].twinx()]
            else:
                axes += [ax[-1]]

        for colors, ax in [(colors_left, axes[0]), (colors_right, axes[-1])]:
            if colors is None:
                continue
            if type(colors) is list:
                c = colors.copy()
                colors = {}
                for m in c:
                    print(str(m))
                    if quantified:
                        if type(m) is str:
                            mol = Molecule(m, verbose=False)
                        else:
                            mol = m
                        color = standard_colors[mol.primary]
                        colors[mol] = color
                    else:
                        color = standard_colors[m]
                        colors[m] = color
            for (key, color) in colors.items():
                if quantified:
                    x, y = get_flux(
                        CV_and_MS,
                        mol=key,
                        tspan=tspan,
                        removebackground=removebackground,
                        background=background,
                        endpoints=endpoints,
                        t_bg=t_bg,
                        unit=unit,
                        verbose=verbose,
                    )
                    if type(key) is not str:
                        key = str(key)  # in case key had been a Molecule object
                    Y_str = key + "_" + unit
                else:

                    Y_str = key + "_" + unit  #
                    x, y = get_signal(
                        CV_and_MS,
                        mass=key,
                        tspan=tspan,
                        removebackground=removebackground,
                        background=background,
                        endpoints=endpoints,
                        t_bg=t_bg,
                        unit=unit,
                        verbose=verbose,
                    )

                try:
                    y_plot = np.interp(
                        t_plot, x, y
                    )  # obs! np.interp has a has a different argument order than Matlab's interp1
                except ValueError:
                    print("x " + str(x) + "\ny " + str(y) + "\nt " + str(t))
                CV_and_MS[
                    Y_str
                ] = y_plot  # add the interpolated value to the dictionary for future use
                # 17C01: but not outside of this function.
                ms_spec = spec.copy()
                if "color" not in ms_spec.keys():
                    ms_spec["color"] = color
                ax.plot(V_plot, y_plot, label=Y_str, **ms_spec)
            if quantified:
                M_str = "cal. signal / [" + unit + "]"
            else:
                M_str = "MS signal / [" + unit + "]"
            ax.set_xlabel(V_str)
            ax.xaxis.set_label_position("top")
            ax.xaxis.tick_top()
            ax.set_ylabel(M_str)
            if leg:
                ax1.legend()

    # if colors_right is not None:
    #    ax[0].set_xlim(ax[1].get_xlim())

    if verbose:
        print("\nfunction 'plot_vs_potential' finished!\n\n")

    for ax in axes:
        if ax is not None:
            ax.tick_params(axis="both", direction="in")  # 17K28

        # parameter order of np.interp is different than Matlab's interp1
    return axes


def plot_vs_potential_EC(
    data,
    tspan=None,
    RE_vs_RHE=None,
    A_el=None,
    cycles="all",
    ax="new",  # spec='k-',
    verbose=True,
    J_str=None,
    V_str=None,
    fig=None,
    spec={},
    t_str=None,
):

    if ax is None or ax == "new":
        fig, ax = plt.subplots()

    V_str, J_str = sync_metadata(
        data, RE_vs_RHE=RE_vs_RHE, A_el=A_el, V_str=V_str, J_str=J_str
    )

    V = data[V_str]
    J = data[J_str]
    if t_str is None:
        if "t_str" in data:
            t_str = data["t_str"]
        else:
            t_str = "time/s"

    # get time variable and plotting indexes
    t = data[t_str]
    if tspan is None:  # then use the whole range of overlap
        try:
            tspan = data["tspan"]
        except KeyError:
            tspan = [min(t), max(t)]

    mask = np.logical_and(tspan[0] < t, t < tspan[-1])

    # plot EC-lab data
    ec_spec = spec.copy()
    print(ec_spec)  # debugging
    if "color" not in ec_spec.keys():
        ec_spec["color"] = "k"

    V_plot, J_plot = V[mask], J[mask]
    ax.plot(V_plot, J_plot, **ec_spec)
    # maybe I should use EC.plot_cycles to have different cycles be different colors. Or rewrite that code here.
    ax.set_xlabel(V_str)
    ax.set_ylabel(J_str)

    if verbose:
        print("\nfunction 'plot_vs_potential_EC' finished!\n\n")
    return ax


def plot_vs_time(Dataset, cols_1="input", cols_2="input", verbose=1):
    """
    Superceded by the more convenient plot_masses and plot_masses_and_I
    """
    if verbose:
        print("\n\nfunction 'plot_vs_time' at your service!")

    if cols_1 == "input":
        data_cols = Dataset["data_cols"]
        prompt = (
            "Choose combinations of time and non-time variables for axis 1, \n"
            + "with every other choice a time variable."
        )
        I_axis_1 = indeces_from_input(data_cols, prompt)
        cols_1 = [
            [data_cols[i], data_cols[j]] for i, j in zip(I_axis_1[::2], I_axis_1[1::2])
        ]

    figure1 = plt.figure()
    axes_1 = figure1.add_subplot(211)
    for pltpair in cols_1:
        label_object = pltpair[1][0:-2]
        if label_object:
            label_string = label_object.group()[:-1]
        else:
            label_string = pltpair[1]
        x = Dataset[pltpair[0]]
        y = np.log(Dataset[pltpair[1]]) / np.log(10)
        axes_1.plot(x, y, label=label_string)

    axes_1.set_xlabel("time / [s]")
    axes_1.set_ylabel("log(MS signal / [a.u.])")
    axes_1.legend()

    if cols_2 == "input":

        data_cols = Dataset["data_cols"]
        prompt = (
            "Choose combinations of time and non-time variables for axis 2, \n"
            + "with every other choice a time variable."
        )
        I_axis_2 = indeces_from_input(data_cols, prompt)
        cols_2 = [
            [data_cols[i], data_cols[j]] for i, j in zip(I_axis_2[::2], I_axis_2[1::2])
        ]

    axes_2 = figure1.add_subplot(212)
    for pltpair in cols_2:
        label_string = pltpair[1]
        x = np.insert(Dataset[pltpair[0]], 0, 0)
        y = np.insert(Dataset[pltpair[1]], 0, 0)
        axes_2.plot(x, y, "k--", label=label_string)
    axes_2.set_ylabel("current / [mA]")
    axes_2.set_xlabel("time / [s]")
    axes_2.legend()
    # so capacitance doesn't blow it up:
    I_plt_top = np.where(x > 2)[0][0]
    y_max = np.max(y[I_plt_top:])
    axes_2.set_ylim(np.min(y), y_max)
    if verbose:
        print("function 'plot_vs_time' finished!\n\n")


def indeces_from_input(options, prompt):
    """something I used all the time back in the (Matlab) days.
        not sure I'll ever actually use it again though"""
    print(
        prompt
        + "\n... enter the indeces you're interested in, in order,"
        + "seperated by spaces, for example:\n>>>1 4 3"
    )
    for nc, option in enumerate(options):
        print(str(nc) + "\t\t " + options[nc])
    choice_string = input("\n")
    choices = choice_string.split(" ")
    choices = [int(choice) for choice in choices]
    return choices


def smooth_data(data, points=3, cols=None, verbose=True):
    """
    Does a moving-average smoothing of data. I don't like it, but
    experencing problems 17G26
    cols should be a list of columns to smooth.
    Operates on the original data set!
    """
    if cols is None:
        cols = data["data_cols"]
    for col in cols:
        if verbose:
            print(
                "smoothening '"
                + col
                + "' with a "
                + str(points)
                + "-point moving average"
            )
        x = data[col].copy()  # in case it's linked to another column.
        c = np.array([1] * points) / points
        # print(str(len(c))) # debugging
        X = np.convolve(x, c, mode="same")
        # convolve doesn't get the endpoints quite right. I can fix them
        for n in range(points):  # fixing endpooints
            X[n] = np.mean(x[0 : n + 1])
            if n > 0:
                X[-n] = np.mean(x[-n:])
        data[col] = X
        # print('len = ' + str(len(x))) # debugging

    return data


def plot_signal(
    MS_data,
    spec={},
    masses="all",
    tspan=None,
    ax="new",
    unit="nA",
    removebackground=None,
    background=None,
    t_bg=None,
    logplot=True,
    saveit=False,
    leg=False,
    name=None,
    override=False,
    verbose=True,
):
    """
    plots selected masses for a selected time range from MS data or EC_MS data
    """
    if name is None:
        try:
            name = MS_data["name"]
        except KeyError:
            try:
                name = MS_data["title"]
            except KeyError:
                name = ""
    if verbose:
        print("\n\nfunction 'plot_signal' at your service! \n Plotting from: " + name)

    if ax == "new":
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    lines = {}

    if masses == "all":  # this is now the default!
        masses = [
            key[:-2] for key in MS_data.keys() if key[0] == "M" and key[-2:] == "-y"
        ]
    elif type(masses) is str:
        masses = [masses]
    if type(masses) is list:
        c = masses
        masses = {}
        for m in c:
            try:
                color = standard_colors[m]
            except KeyError:
                print("Waring: no standard color for " + m + ". Using black.")
                color = "k"
            masses[m] = color

    for mass, color in masses.items():
        if verbose:
            print("plotting: " + mass)
        try:
            x, y = get_signal(
                MS_data,
                mass,
                unit=unit,
                tspan=tspan,
                override=override,
                verbose=verbose,
                removebackground=removebackground,
                t_bg=t_bg,
                background=background,
            )
            if len(x) == 0:
                print("WARNING: no data for " + mass)
                continue
        except KeyError:
            print("WARNING: Can't get signal for " + str(mass))
            continue
        if len(x) == 0:
            print(
                "WARNING: get_signal returned vector of zero length for "
                + mass
                + ". plot_signal is skipping that mass."
            )
            continue
        lines[mass] = ax.plot(x, y, color, label=mass, **spec)
        # as it is, lines is not actually used for anything
    if leg:
        if type(leg) is not str:
            leg = "lower right"
        ax.legend(loc=leg)
    ax.set_xlabel("time / [s]")
    ax.set_ylabel("MS signal / [" + unit + "]")
    if logplot:
        ax.set_yscale("log")
    ax.tick_params(axis="both", direction="in")  # 17K28
    if verbose:
        print("function 'plot_signal' finsihed! \n\n")
    return ax


def plot_masses(*args, **kwargs):
    print("plot_masses renamed plot_signal. Remember that next time!")
    return plot_signal(*args, **kwargs)


def plot_flux(
    MS_data,
    mols={"H2": "b", "CH4": "r", "C2H4": "g", "O2": "k"},
    tspan=None,
    ax="new",
    removebackground=False,
    background="constant",
    endpoints=5,
    t_bg=None,
    A_el=None,
    unit="nmol/s",
    smooth_points=0,
    logplot=True,
    leg=False,
    spec={},
    override=False,
    verbose=True,
):
    """
    Plots the molecular flux to QMS in nmol/s for each of the molecules in
    'fluxes.keys()', using the primary mass and the F_cal value read from
    the molecule's text file, with the plot specs from 'fluxes.values()'
    """
    if verbose:
        print("\n\nfunction 'plot_flux' at your service!\n")
    if ax == "new":
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    # note, tspan is processed in get_flux, and not here!

    if type(mols) is not list and type(mols) is not dict:
        # then it's probably just one molecule object.
        mols = [mols]

    if type(mols) is list:
        c = mols
        mols = {}
        for m in c:
            if type(m) is str:
                mol = Molecule(m, verbose=False)
            else:
                mol = m  # this function should accept a list of Molecule instances!
            color = standard_colors[mol.primary]
            mols[mol] = color
            print("mol={}, primary={}, color={}".format(mol.name, mol.primary, color))

    for (mol, color) in mols.items():
        try:
            [x, y] = get_flux(
                MS_data,
                mol,
                unit=unit,
                tspan=tspan,
                removebackground=removebackground,
                background=background,
                t_bg=t_bg,
                endpoints=endpoints,
                override=override,
                verbose=verbose,
            )
            if smooth_points:
                y = smooth(y, smooth_points)
        except KeyError:
            print("Can't get signal for " + str(mol))
            continue
        if type(mol) is str:
            l = mol
        else:
            l = mol.name
        print("color={}".format(color))  # debugging
        ax.plot(x, y, color=color, label=l, **spec)
    if leg:
        if type(leg) is not str:
            leg = "lower right"
        ax.legend(loc=leg)
    ax.set_xlabel("time / [s]")
    ylabel = "cal. signal / [" + unit + "]"

    ax.set_ylabel(ylabel)
    if logplot:
        ax.set_yscale("log")

    ax.tick_params(axis="both", direction="in")  # 17K28

    if verbose:
        print("\nfunction 'plot_flux' finished!\n\n")
    return ax


def plot_experiment_EC(
    data,
    tspan=None,
    verbose=True,
    RE_vs_RHE=None,
    A_el=None,
    ax="new",
    # mols will overide masses will overide colors
    V_color=None,
    J_color=None,
    V_label=None,
    J_label=None,
    t_str=None,
    J_str=None,
    V_str=None,
    fig=None,
    spec={},
):
    if ax == "new":
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax = [ax1, ax2]

    # --------- get tspan, V_str, and J_str from input and/or dataset -------- #
    if tspan is None:  # then use the range of overlap
        try:
            tspan = data["tspan"]  # changed from 'tspan_2' 17H09
        except KeyError:
            tspan = "all"
    if type(tspan) is str and not tspan == "all":
        tspan = data[tspan]

    if t_str is None:
        if "t_str" in data:
            t_str = data["t_str"]
        else:
            t_str = "time/s"
    if V_str is None or J_str is None or RE_vs_RHE is not None or A_el is not None:
        V_str_0, J_str_0 = sync_metadata(
            data, RE_vs_RHE=RE_vs_RHE, A_el=A_el, verbose=verbose
        )
        # added 16J27... problem caught 17G26, fixed in sync_metadata
    if V_str is None:  # this way I can get it to plot something other than V and J.
        V_str = V_str_0
    if J_str is None:
        J_str = J_str_0

    if A_el in data:
        A_el = data["A_el"]
    else:
        A_el = 1

    # ---------- make sure I can plot the electrochemistyr data --------- #
    plotpotential = True
    plotcurrent = True
    try:
        t = data[t_str]
    except KeyError:
        print(
            "data doesn't contain '" + str(t_str) + "', i.e. t_str. Can't plot EC data."
        )
        plotpotential = False
        plotcurrent = False
    try:
        V = data[V_str]
    except KeyError:
        print(
            "data doesn't contain '"
            + str(V_str)
            + "', i.e. V_str. Can't plot that data."
        )
        plotpotential = False
    try:
        J = data[J_str]
    except KeyError:
        print(
            "data doesn't contain '"
            + str(J_str)
            + "', i.e. J_str. Can't plot that data."
        )
        plotcurrent = False

    # -------- cut the electrochemistry data according to tspan ------ #
    if type(tspan) is not str and (plotcurrent or plotpotential):
        mask = np.logical_and(tspan[0] < t, t < tspan[-1])
        t = t[mask]
        # print(np.where(mask)) #debugging
        if plotpotential:
            V = V[mask]
        if plotcurrent:
            J = J[mask]

    # ---------- and plot the electrochemistry data! --------------- #
    if plotcurrent:
        if plotpotential:
            i_ax = 1
        else:
            i_ax = 0
        ax[i_ax].plot(t, J, color=J_color, label=J_label, **spec)
        ax[i_ax].set_ylabel(J_str)
        ax[i_ax].set_xlabel("time / [s]")
        xlim = ax[i_ax - 1].get_xlim()
        ax[i_ax].set_xlim(xlim)

        if i_ax == 2:
            colorax(ax[i_ax], J_color, "right")
        else:
            colorax(ax[i_ax], J_color, "left")
        ax[i_ax].tick_params(axis="both", direction="in")  # 17K28

    if plotpotential:
        i_ax = 0
        ax[i_ax].plot(t, V, color=V_color, label=V_label, **spec)
        ax[i_ax].set_ylabel(V_str)
        xlim = ax[i_ax - 1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        colorax(ax[i_ax], V_color, "left")
        ax[i_ax].tick_params(axis="both", direction="in")  # 17K28

    ax[0].set_xlabel(t_str)
    if plotcurrent or plotpotential:
        ax[0].xaxis.set_label_position("top")
        ax[1].set_xlabel(t_str)
        if tspan is not None and not type(tspan) is str:
            ax[1].set_xlim(tspan)

    # -------- finishing up -------- #
    print("\nfunction 'plot_experiment_EC' finished!\n\n")
    return ax


def plot_experiment(
    EC_and_MS,
    colors=None,
    tspan=None,
    tspan_EC=None,
    overlay=False,
    logplot=[True, False],
    verbose=True,
    plotpotential=True,
    plotcurrent=True,
    ax="new",
    emphasis="ms",
    RE_vs_RHE=None,
    A_el=None,
    removebackground=None,
    background=None,
    endpoints=5,
    t_bg=None,
    saveit=False,
    name=None,
    leg=False,
    unit=None,
    masses="all",
    masses_left=None,
    masses_right=None,
    mols=None,
    mols_left=None,
    mols_right=None,
    # mols will overide masses will overide colors
    V_color="k",
    J_color="0.5",
    V_label=None,
    J_label=None,
    t_str=None,
    J_str=None,
    V_str=None,
    fig=None,
    return_fig=False,
    smooth_points=0,
    override=False,
    spec={},
    hspace=0.1,
    left_space=0.15,
    right_space=0.85,  # for gridspec
):
    """
    TODO: write proper documentation!
    this plots signals or fluxes on one axis and current and potential on other axesaxis

    background can be:
        'constant' - subtracts minimum value or given value
        'linear' - linear background extrapolated between endpoints
        'preset' - finds backgrounds in molecule objects.
    By default, no background is subtracted.
    """
    if name is None:
        try:
            name = EC_and_MS["name"]
        except KeyError:
            try:
                name = EC_and_MS["title"]
            except KeyError:
                name = ""

    if verbose:
        print(
            "\n\nfunction 'plot_experiment' at your service!\n Plotting from: " + name
        )

    if "data_type" in EC_and_MS and EC_and_MS["data_type"][0:2] == "EC":
        ax = plot_experiment_EC(
            EC_and_MS,
            tspan=tspan,
            verbose=verbose,
            RE_vs_RHE=RE_vs_RHE,
            A_el=A_el,
            ax=ax,
            # mols will overide masses will overide colors
            V_color=V_color,
            J_color=J_color,
            V_label=V_label,
            J_label=J_label,
            t_str=t_str,
            J_str=J_str,
            V_str=V_str,
            fig=None,
            spec={},
        )
        return ax

    # ----------- prepare the axes on which to plot ------------ #
    if ax == "new":
        if fig is None:
            figure1 = plt.figure()
        else:
            figure1 = fig
            plt.figure(figure1.number)
            print("plot_expeiriment using " + str(fig))
        if overlay:
            ax = [figure1.add_subplot(111)]
            ax += [ax[0].twinx()]
        else:
            if emphasis == "MS":
                gs = gridspec.GridSpec(12, 1)
                # gs.update(hspace=0.025)
                ax = [plt.subplot(gs[0:8, 0])]
                ax += [plt.subplot(gs[8:12, 0])]
            elif emphasis == "ms":
                gs = gridspec.GridSpec(5, 1)
                # gs.update(hspace=0.025)
                ax = [plt.subplot(gs[0:3, 0])]
                ax += [plt.subplot(gs[3:5, 0])]
            elif emphasis == "EC":
                gs = gridspec.GridSpec(3, 1)
                # gs.update(hspace=0.025)
                ax = [plt.subplot(gs[0, 0])]
                ax += [plt.subplot(gs[1:3, 0])]
            else:
                gs = gridspec.GridSpec(8, 1)
                # gs.update(hspace=0.025)
                ax = [plt.subplot(gs[0:4, 0])]
                ax += [plt.subplot(gs[4:8, 0])]
            if plotcurrent and plotpotential:
                ax += [ax[1].twinx()]
                ax[1].set_zorder(ax[2].get_zorder() + 1)  # doesn't work
                ax[1].patch.set_visible(False)  # hide the 'canvas'
            gs.update(hspace=hspace, left=left_space, right=right_space)

    # --------- get tspan, V_str, and J_str from input and/or dataset -------- #
    if tspan is None:  # then use the range of overlap
        try:
            tspan = EC_and_MS["tspan"]  # changed from 'tspan_2' 17H09
        except KeyError:
            tspan = "all"
    if type(tspan) is str and not tspan == "all":
        tspan = EC_and_MS[tspan]
    if type(logplot) is not list:
        logplot = [logplot, False]

    if t_str is None:
        if "t_str" in EC_and_MS:
            t_str = EC_and_MS["t_str"]
        else:
            t_str = "time/s"
    if V_str is None or J_str is None or RE_vs_RHE is not None or A_el is not None:
        V_str_0, J_str_0 = sync_metadata(
            EC_and_MS, RE_vs_RHE=RE_vs_RHE, A_el=A_el, verbose=verbose
        )
        # added 16J27... problem caught 17G26, fixed in sync_metadata
    if V_str is None:  # this way I can get it to plot something other than V and J.
        V_str = V_str_0
    if J_str is None:
        J_str = J_str_0

    if A_el in EC_and_MS:
        A_el = EC_and_MS["A_el"]
    else:
        A_el = 1

    # ----------- parse input on which masses / fluxes to plot ------- #

    #    print(masses)
    quantified = False
    # print(type(colors))
    # if type(colors) is list and type(colors[0]) is not str:
    #    print(type(colors[0]))
    if mols is not None:
        quantified = True
    elif (
        (type(colors) is dict and list(colors.keys())[0][0] == "M")
        or (type(colors) is list and type(colors[0]) is str and colors[0][0] == "M")
        or (type(colors) is str and colors[0] == "M")
        or colors is None
    ):
        if verbose:
            print("uncalibrated data to be plotted.")
        if masses is None:
            masses = colors
    else:
        quantified = True
        mols = colors

    if not quantified and masses == "all":  # this is now the default!
        masses = [
            key[:-2] for key in EC_and_MS.keys() if key[0] == "M" and key[-2:] == "-y"
        ]
    print("quantified = " + str(quantified))  # debugging

    if removebackground is None and background is None and t_bg is None:
        removebackground = False

    if type(mols) is dict:
        mols = list(mols.values())
    try:
        if type(mols[0]) in (list, tuple):
            mols_left = mols[0]
            mols_right = mols[1]
            mols = mols_left
    except (IndexError, TypeError):
        pass
    if mols_left is not None and mols is None:
        mols = mols_left
    try:
        if type(masses[0]) in (list, tuple):
            masses_left = masses[0]
            masses_right = masses[1]
            masses = masses_left
    except (IndexError, TypeError):
        pass
    if masses_left is not None and masses is None:
        masses = masses_left
    if removebackground == "right":
        removebackground_right = True
        removebackground_left = False
    elif removebackground == "left":
        removebackground_left = True
        removebackground_right = False
    else:
        removebackground_left = removebackground
        removebackground_right = removebackground

    # ----------- Plot the MS signals! ------------- #
    if quantified:
        if unit is None:
            unit = "pmol/s"
        print("removebackground = " + str(removebackground))  # debugging
        plot_flux(
            EC_and_MS,
            mols=mols,
            tspan=tspan,
            A_el=A_el,
            spec=spec,
            ax=ax[0],
            leg=leg,
            logplot=logplot[0],
            unit=unit,
            removebackground=removebackground_left,
            background=background,
            endpoints=endpoints,
            t_bg=t_bg,
            override=override,
            smooth_points=smooth_points,
            verbose=verbose,
        )
        if mols_right is not None:
            ax += [ax[0].twinx()]
            plot_flux(
                EC_and_MS,
                mols=mols_right,
                tspan=tspan,
                A_el=A_el,
                spec=spec,
                ax=ax[-1],
                leg=leg,
                logplot=logplot[0],
                unit=unit,
                removebackground=removebackground_right,
                background=background,
                endpoints=endpoints,
                t_bg=t_bg,
                override=override,
                smooth_points=smooth_points,
                verbose=verbose,
            )
    else:
        if unit is None:
            unit = "pA"
        plot_signal(
            EC_and_MS,
            masses=masses,
            tspan=tspan,
            spec=spec,
            ax=ax[0],
            leg=leg,
            logplot=logplot[0],
            unit=unit,
            override=override,
            verbose=verbose,
            removebackground=removebackground_left,
            background=background,
            t_bg=t_bg,
        )
        if masses_right is not None:
            ax += [ax[0].twinx()]
            plot_signal(
                EC_and_MS,
                masses=masses_right,
                tspan=tspan,
                spec=spec,
                ax=ax[-1],
                leg=leg,
                logplot=logplot[0],
                unit=unit,
                override=override,
                verbose=verbose,
                removebackground=removebackground_right,
                background=background,
                t_bg=t_bg,
            )
    if not overlay:
        ax[0].set_xlabel("")
        ax[0].xaxis.tick_top()

    if tspan is not None and not type(tspan) is str:
        ax[0].set_xlim(tspan)

    # ---------- make sure I can plot the electrochemistyr data --------- #

    try:
        t = EC_and_MS[t_str]
    except KeyError:
        print(
            "data doesn't contain '" + str(t_str) + "', i.e. t_str. Can't plot EC data."
        )
        plotpotential = False
        plotcurrent = False
    try:
        V = EC_and_MS[V_str]
    except KeyError:
        print(
            "data doesn't contain '"
            + str(V_str)
            + "', i.e. V_str. Can't plot that data."
        )
        plotpotential = False
    try:
        J = EC_and_MS[J_str]
    except KeyError:
        print(
            "data doesn't contain '"
            + str(J_str)
            + "', i.e. J_str. Can't plot that data."
        )
        plotcurrent = False

        # to check if I have problems in my dataset
    #    print('len(t) = ' + str(len(t)) +
    #          '\nlen(V) = ' + str(len(V)) +
    #          '\nlen(J) = ' + str(len(J)))
    # print(tspan) # debugging

    # -------- cut the electrochemistry data according to tspan ------ #
    if tspan_EC is None:
        tspan_EC = tspan
    if type(tspan_EC) is not str and (plotcurrent or plotpotential):
        mask = np.logical_and(tspan_EC[0] < t, t < tspan_EC[-1])
        t = t[mask]
        # print(np.where(mask)) #debugging
        if plotpotential:
            V = V[mask]
        if plotcurrent:
            J = J[mask]

    # ---------- and plot the electrochemistry data! --------------- #
    if plotcurrent:
        if plotpotential:
            i_ax = 2
        else:
            i_ax = 1
        ax[i_ax].plot(t, J, color=J_color, label=J_label, **spec)
        ax[i_ax].set_ylabel(J_str)
        ax[i_ax].set_xlabel("time / [s]")
        xlim = ax[i_ax - 1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        if logplot[1]:
            ax[i_ax].set_yscale("log")

        if i_ax == 2:
            colorax(ax[i_ax], J_color, "right")
        else:
            colorax(ax[i_ax], J_color, "left")
        ax[i_ax].tick_params(axis="both", direction="in")  # 17K28

    if plotpotential:
        i_ax = 1
        ax[i_ax].plot(t, V, color=V_color, label=V_label, **spec)
        ax[i_ax].set_ylabel(V_str)
        if len(logplot) > 2:
            if logplot[2]:
                ax[i_ax].set_yscale("log")
        xlim = ax[i_ax - 1].get_xlim()
        ax[i_ax].set_xlim(xlim)
        colorax(ax[i_ax], V_color, "left")
        ax[i_ax].tick_params(axis="both", direction="in")  # 17K28

    ax[0].set_xlabel(t_str)
    if plotcurrent or plotpotential:
        ax[0].xaxis.set_label_position("top")
        ax[1].set_xlabel(t_str)
        if tspan is not None and not type(tspan) is str:
            ax[1].set_xlim(tspan)

    # -------- finishing up -------- #
    if saveit:
        figure1.savefig(name + ".png")

    # if colors_right is not None:
    #    ax[0].set_xlim(ax[1].get_xlim())   # probably not necessary

    if return_fig and (fig is None):
        fig = ax[0].get_figure()
    if verbose:
        print("function 'plot_experiment' finished!\n\n")
    if return_fig:
        return fig, ax
    else:
        return ax


def plot_masses_and_I(*args, **kwargs):
    print(
        "\n\n'plot_masses_and_I' has been renamed 'plot_experiment'. Remember that next time!"
    )
    return plot_experiment(*args, **kwargs)


def plot_folder(
    folder_name,
    colors={"M2": "b", "M4": "r", "M18": "0.5", "M28": "g", "M32": "k"},
    RE_vs_RHE=None,
    A_el=None,
):
    """
    Plots an EC and MS data from an entire folder, generally corresponding to
    a full day of measurements on one sample.
    Will probably only be used to get an overview.
    Could add text showing starts of the data files
    """
    from .Data_Importing import import_folder
    from .Combining import synchronize

    Datasets = import_folder(folder_name)
    Combined_data = synchronize(Datasets, t_zero="first", append=True)
    sync_metadata(Combined_data, RE_vs_RHE, A_el)
    return plot_experiment(Combined_data, colors=colors)


def plot_lines_x(values, ax="new", ylims=None, **kwargs):
    if ax == "new":
        fig, ax = plt.subplots()
    if ylims is None:
        ylims = ax.get_ylim()
        ax.set_ylim(ylims)
    for v in values:
        ax.plot([v, v], ylims, **kwargs)
    return ax


def plot_datapoints(
    integrals,
    colors,
    ax="new",
    label="",
    X=None,
    X_str="V",
    logplot=True,
    specs={},
    Xrange=None,
):
    """
    integrals will most often come from functino 'get_datapoitns' in module
    Integrate_Signals
    """
    if ax == "new":
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    if X is None:
        X = integrals[X_str]

    for (quantity, color) in colors.items():
        # Here I just assme they've organized the stuff right to start with.
        # I could alternately use the more intricate checks demonstrated in
        # DataPoints.plot_errorbars_y
        value = integrals[quantity]
        if type(Xrange) is dict:
            Xrange_val = Xrange[quantity]
        else:
            Xrange_val = Xrange
        if type(color) is dict:
            plot_datapoints(
                value,
                color,
                ax=ax,
                logplot=logplot,
                label=label + quantity + "_",
                X=X,
                Xrange=Xrange_val,
                specs=specs,
            )
        else:
            if type(color) is tuple:  # note a list can be a color in rbg
                spec = color[0]
                color = color[1]
                if "markersize" not in specs:
                    specs["markersize"] = 5
            else:
                spec = "."
                if "markersize" not in specs:
                    specs["markersize"] = 15
            # print(quantity + '\n\tvalue=' + str(value) +
            #        '\n\tcolor=' + str(color) + '\n\tV=' + str(V))
            # print(quantity + ' ' + str(color))
            if Xrange is not None:
                I_keep = np.array(
                    [
                        I
                        for (I, X_I) in enumerate(X)
                        if Xrange_val[0] <= float(np.round(X_I, 2)) <= Xrange_val[1]
                    ]
                )
                X_plot = np.array(X)[I_keep]
                value_plot = np.array(value)[I_keep]
                # there was a mindnumbing case of linking here.
                # tried fix it with .copy(), but new variable names needed.
            else:
                X_plot = X
                value_plot = value
            ax.plot(
                X_plot, value_plot, spec, color=color, label=label + quantity, **specs,
            )
    if logplot:
        ax.set_yscale("log")
    return ax


def plot_operation(
    cc=None,
    t=None,
    j=None,
    z=None,
    tspan=None,
    results=None,
    plot_type="heat",
    ax="new",
    colormap="inferno",
    aspect="auto",
    unit="pmol/s",
    color="g",
    colorbar=True,
    dimensions=None,
    verbose=True,
):
    if verbose:
        print("\n\nfunction 'plot_operation' at your service!\n")
    # and plot!
    if type(cc) is dict and results is None:
        results = cc
        cc = None
    if results is None:
        results = {}  # just so I don't get an error later
    if cc is None:
        cc = results["cc"]
    if t is None:
        if "t" in results:
            t = results["t"]
        elif "x" in results:
            t = results["x"]
        else:
            t = np.linspace(0, 1, np.size(cc, axis=0))
    if z is None:
        if "z" in results:
            z = results["z"]
        elif "y" in results:
            z = results["y"]
        else:
            z = np.linspace(0, 1, np.size(cc, axis=1))
    if j is None:
        if j in results:
            j = results["j"]
        else:
            j = cc[0, :]
    if dimensions is None:
        if "dimensions" in results:
            dimensions = results["dimensions"]
        else:
            dimensions = "tz"

    if tspan is None:
        tspan = [t[0], t[-1]]
    if plot_type == "flux":
        if ax == "new":
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
        else:
            ax1 = ax  # making a heat map will only work with a new axis.
        ax1.plot(t, j, label="simulated flux")
        ax1.set_xlabel("time / [s]")
        ax1.set_ylabel("flux / [" + unit + "]")
        axes = ax1

    elif plot_type == "heat" or plot_type == "both":
        if ax == "new":
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
        elif type(ax) is list:
            ax1 = ax[0]
        else:
            ax1 = ax

        # t_mesh, x_mesh = np.meshgrid(t,x)
        # img = ax1.contourf(t_mesh, x_mesh*1e6, np.transpose(cc,[1,0]), cmap='Spectral',
        #                   levels=np.linspace(np.min(cc),np.max(cc),100))

        # imshow objects seem more versatile than contourf for some reason.

        trange = [min(t), max(t)]
        if dimensions[0] == "x":
            trange = [t * 1e3 for t in trange]  # m to mm
        zrange = [min(z * 1e6), max(z * 1e6)]

        img = ax1.imshow(
            np.transpose(cc, [1, 0]),
            extent=trange[:] + zrange[:],  # have to be lists here!
            aspect=aspect,
            origin="lower",
            cmap=colormap,
        )

        #        divider = make_axes_locatable(ax1)
        #        cax = divider.append_axes("right", size="5%", pad=0.05)
        # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph

        if colorbar:
            cbar = plt.colorbar(img, ax=ax1)
            cbar.set_label("concentration / [mM]")
            if dimensions[0] == "t":
                ax1.set_xlabel("time / [s]")
            elif dimensions[0] == "x":
                ax1.set_xlabel("position / [mm]")
            ax1.set_ylabel("position / [um]")

        #        print('plot_type = ' + plot_type)
        if plot_type == "both":
            if type(ax) is list:
                ax2 = ax[1]
            else:
                ax2 = ax1.twinx()
            ax2.set_ylabel("flux / [" + unit + "]")
            ax2.plot(t, j, "-", color=color)
            cbar.remove()
            ax3 = img.figure.add_axes([0.85, 0.1, 0.03, 0.8])
            cbar = plt.colorbar(img, cax=ax3)
            cbar.set_label("concentration / [mM]")
            ax1.set_xlim(tspan)
            print("returning three axes!")
            axes = [ax1, ax2, ax3]
        elif colorbar:
            axes = [ax1, cbar]
        else:
            axes = ax1
    if verbose:
        print("\nfunction 'plot_operation' finished!\n\n")
    return axes


def set_figparams(
    figwidth=8, aspect=4 / 3, fontsize=7, figpad=0.15, figgap=0.08, style=None
):
    import matplotlib as mpl

    if style is not None:
        mpl.style.use(style)

    # figwidth=8  #figwidth in cm width 20.32cm = 8inches being standard and thesis textwidth being 12.
    # aspect=4/3  #standard is 4/3

    # fontsize=7  #standard is 12.0pt, thesis is 10.0pt and footnotesize is 8.0pt and lower case seems to be 2/3 of full font size, which makes 7pt "nice" for thesis plotting

    realfigwidth = (
        20 * (fontsize / 12) * 1.2
    )  # a factor 1.2 makes all lines "equaly thick" - make this 1 for figwidth=4cm (sniff2fig6)
    # figsize=[20.32,15.24]
    figsize = [realfigwidth, realfigwidth / aspect]

    mpl.rc("font", size=fontsize * (realfigwidth / figwidth))

    mpl.rc(
        "mathtext",
        fontset="custom",
        rm="Helvetica",
        # it='Helvetica:italic',
        # bf='Helvetica:bold',
    )

    mpl.rc(
        "figure",
        figsize=[figsize[0] / 2.54, figsize[1] / 2.54],
        dpi=100 * 2.54 * figwidth / realfigwidth,
    )

    # figpad=0.14  #fraction of figure size
    # figgap=0.08  #fraction of figure size
    mpl.rc(
        "figure.subplot",
        left=figpad,
        right=1 - figpad,
        bottom=figpad,
        top=1 - figpad,
        hspace=figgap,
    )
    mpl.rc("xtick", labelsize="small")
    mpl.rc("ytick", labelsize="small")

    # mpl.rc('axes', labelweight='medium')

    mpl.rc("savefig", dpi=250 * 2.54 * figwidth / realfigwidth)
