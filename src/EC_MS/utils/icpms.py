#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:33:05 2019

@author: scott
"""
import re
import numpy as np
from matplotlib import pyplot as plt


unit_to_ppb = {
    "ppt": 1e-3,
    "ng/l": 1e-3,
    "ppb": 1,
    "ug/l": 1,
    "ppm": 1000,
    "mg/l": 1000,
}
float_match = "[-]?\d+[\.]?\d*(e[-]?\d+)?"  # matches floats like '-3.5e4' or '7' or '245.13' or '1e-15'


def get_calibration_function(
    ICPMS_data,
    mass,
    unit="ug/l",
    units=[],
    ax=None,
    wash_str=None,
    plot_unit=None,
    factor=1,
):
    """
    Given ICPMS_data, a dictionary of the structure returned by load_ICPMS_data,
    this function finds the calibration standards (which have the specified unit,
    'ug/l' by default, in their name), parses their name to get the concentration,
    and plots a calibration curve. It then returns the function given in that curve.
    """
    if len(units) == 0:
        units = [unit]

    calibration_keys = [
        key for key in ICPMS_data.keys() if len([u for u in units if u in key]) > 0
    ]

    # print(calibration_keys) # debugging

    cal_amounts = np.array([])
    cal_counts = np.array([])

    for key in calibration_keys:
        unit_i = next(u for u in units if u in key)
        try:
            number = re.search(float_match, key).group()
            amount_i = float(number)
        except (AttributeError, ValueError):
            print("WARNING: could't match a float in '" + key + "'. skipping.")
            continue
        amount = amount_i * unit_to_ppb[unit_i] / unit_to_ppb[unit]
        cal_amounts = np.append(cal_amounts, amount)

        cal_counts = np.append(cal_counts, ICPMS_data[key][mass])
        # print(str(cal_amounts) + '\n' + str(cal_counts)) # debugging

    ln_cal_amounts = np.log(cal_amounts)
    ln_cal_counts = np.log(cal_counts)

    p = np.polyfit(ln_cal_counts, ln_cal_amounts, deg=1)
    print(p)  # debugging

    def counts_to_amount(counts):
        ln_counts = np.log(counts)
        ln_ppb = p[0] * ln_counts + p[1]
        ppb = np.exp(ln_ppb)
        return ppb

    if ax == "new":
        fig, ax = plt.subplots()
    if ax is not None:
        ax.plot(cal_amounts, cal_counts, "ks", markersize=7)
        x_fit, y_fit = counts_to_amount(cal_counts), cal_counts

        if wash_str is not None:
            wash = ICPMS_data[wash_str][mass]
            mean = np.mean(wash)
            std = np.std(wash)

            y_fit = np.append(mean, y_fit)
            x_fit = np.append(counts_to_amount(mean), x_fit)

        ax.plot(x_fit, y_fit, "r--")
        if wash_str is not None:
            xlim = ax.get_xlim()
            ax.plot(xlim, [mean, mean], "k--")
            ax.plot(xlim, [mean + 3 * std, mean + 3 * std], "k:")
            # ax.plot(xlim, [mean-3*std, mean-3*std], 'k:')
            ax.set_xlabel("amount / [" + unit + "]")
            ax.set_ylabel("counts")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(xlim)
        if plot_unit is not None:
            ax2 = ax.twiny()
            ax2.set_xlim([l * factor for l in ax.get_xlim()])
            ax2.set_xscale("log")
            ax2.set_xlabel("amount / [" + plot_unit + "]")

    return counts_to_amount


def load_ICPMS_data(ICPMS_file, sep="\t"):
    """
    This assumes a simple export, including only the averaged intensities.
    If it doesn't work, you're probably included extra stuff when exporting
    from the ICP-MS software!
    Results are returned as a nested dictionary with the sample name as the
    outer key and the mass as the inner key. If a sample name is repeated,
    the function appends the values.
    """

    data = {"file": ICPMS_file}
    data_start = False
    header = ""

    with open(ICPMS_file, "r") as f:

        while True:
            line = f.readline().rstrip()
            if len(line) == 0:
                # ^ this is the case when we've reached the end of the file
                break

            if not data_start:
                header += line
                if "(KED)" in line:
                    cs = line.split(sep)
                    masses = {}
                    for i, c in enumerate(cs[2:]):
                        # there are two blank columns at the start
                        mass = c.split(" ")[0]
                        masses[i] = mass
                    data["masses"] = masses
                    # print('cs = ' + str(cs) + '\nmasses = ' + str(masses)) # debugging

                if "Y (cps)" in line:
                    data_start = True
                    continue

            if data_start:
                cs = line.split(sep)
                # the first column is a number that (I think) isn't needed
                try:
                    name = cs[1]
                except IndexError:
                    print("WARNING! error on line : \n" + line)
                if not name in data:
                    data[name] = {}
                for i, c in enumerate(cs[2:]):
                    try:
                        n = float(c)
                    except:
                        print("error on line:\n\t" + line)
                        print("can't convert " + c + " to float.")
                        raise
                    mass = masses[i]
                    if mass in data[name]:
                        data[name][mass] = np.append(data[name][mass], n)
                    else:
                        data[name][mass] = n

    return data
