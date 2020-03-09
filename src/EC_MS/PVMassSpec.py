# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 11:06:17 2020

@author: Spectro Inlets Lab_1
"""

import re
import numpy as np

from .Data_Importing import timestring_to_epoch_time

PVMassSpec_time_match = r"[0-9]{2}-[0-9]{2}-[0-9]{4} [0-9]{2}'[0-9]{2}'[0-9]{2}"  # it has a ridiculous format


def PVMS_title_to_timestring(title):
    match = re.search(PVMassSpec_time_match, title)
    if not match:
        print("Warning\! no timestamp found in " + title)
        return
    match_str = match.group()
    m, d, Y = match_str[0:2], match_str[3:5], match_str[6:10]
    H, M, S = match_str[11:13], match_str[14:16], match_str[17:19]

    timestring = f"{Y}/{m}/{d} {H}:{M}:{S}"
    return timestring


def read_PVMS(file_path, delim="\t", t_str="Time Relative (sec)"):

    with open(file_path, "r") as f:
        lines = f.readlines()

    data = {
        "file": file_path,
        "header": "",
        "timecols": {},
    }

    N_col_head = len(
        lines
    )  # this will decrease when the loop knows when the column header line is comming

    nondata_cols = []  # this will store abstime, so I don't have to parse.
    for n, line in enumerate(lines):
        l = line.strip()
        if n < N_col_head:
            data["header"] = data["header"] + line
            if len(l) == 0:
                N_col_head = n + 1
        elif n == N_col_head:
            data_cols = l.split(delim)
            for col in data_cols:
                data[col] = np.array([])
                data["timecols"][col] = "Time Relative (sec)"
        elif n > N_col_head:
            for col, val in zip(data_cols, l.split(delim)):
                if col in nondata_cols:
                    data[col] += [val]
                    continue
                try:
                    try:
                        x = eval(val)
                    except NameError:
                        if "nan" in val:
                            x = np.nan
                        else:
                            raise SyntaxError
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

    data["data_cols"] = set([col for col in data_cols if col not in nondata_cols])
    data["t_str"] = t_str
    if t_str not in data["data_cols"]:
        print("Warning!!! " + t_str + " not in data from " + file_path)
    else:
        data["timecols"] = dict([(col, t_str) for col in data["data_cols"]])

    timestring = PVMS_title_to_timestring(file_path)
    if not timestring:
        print(
            "couldn't find the timestring in "
            + file_path
            + ". Using the first point in Time Absolute."
        )
    if timestring:
        # print(timestring + '!!!') # debugging
        tstamp = timestring_to_epoch_time(timestring)
    else:
        try:
            timestring = data["Time Absolute (Date_Time)"][0]
        except (KeyError, IndexError):
            print("Warning!!! could't get timestring from Time Absolute either.")
        try:
            tstamp = data["Time Absolute (UTC)"][0]
        except (KeyError, IndexError):
            tstamp = timestring_to_epoch_time(timestring)

    data["timestring"] = timestring
    data["tstamp"] = tstamp

    rename_PVMS_cols(data)

    return data


def read_PVMS_spectrum(*args, index=0, **kwargs):
    """
    calls read_PVMS_spectra with args and kwargs, then returns spectra[index].
    """
    spectra = read_PVMS_spectra(*args, **kwargs)
    spectrum = spectra[index]
    return spectrum


def rename_PVMS_cols(data):
    data_cols = data["data_cols"].copy()
    timecol = "Time Relative (sec)"
    for col in data_cols:
        if "_amu" in col:
            mass = "M" + col.split("_amu")[0]
            xcol = mass + "-x"
            ycol = mass + "-y"
            data[ycol] = data[col]
            data[xcol] = data[timecol]
            data["data_cols"].add(xcol)
            data["data_cols"].add(ycol)
            data["timecols"][ycol] = xcol
            data["timecols"][xcol] = xcol
    # remove_negatives(dataset) # not sure if this is useful.
    return data  # not really necessary
