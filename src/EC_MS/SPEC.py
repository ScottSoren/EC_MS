# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:58:57 2020

@author: scott


This module contains stuff that was copied into Data_Importing from EC_Xray.
It should probably just be removed from EC_MS altogether.
"""


import re

from .parsing_tools import (
    get_creation_timestamp,
    timetag_to_timestamp,
    get_empty_set,
    numerize,
    remove_comments,
    float_match,
)


def load_from_csv(filepath, multiset=False, timestamp=None, verbose=True):
    """
    This function is made a bit more complicated by the fact that some csvs
    seem to have multiple datasets appended, with a new col_header line as the
    only indication. If multiset=True, this will separate them and return them
    as a list.
    if timestamp = None, the timestamp will be the date created
    I hate that SPEC doesn't save absolute time in a useful way.
    """
    if verbose:
        print("function 'load_from_csv' at your service!")
    if timestamp is None:
        a = re.search("[0-9]{2}h[0-9]{2}", filepath)
        if a is None:
            print("trying to read creation time")
            timestamp = get_creation_timestamp(filepath)
        else:
            print("getting timestamp from filename " + filepath)
            timestamp = timetag_to_timestamp(filepath)

    with open(filepath, "r") as f:  # read the file!
        lines = f.readlines()
    colheaders = [col.strip() for col in lines[0].split(",")]
    data = get_empty_set(
        set(colheaders), title=filepath, timestamp=timestamp, data_type="SPEC"
    )
    datasets = []

    for line in lines[1:]:  # put data in lists!
        vals = [val.strip() for val in line.split(",")]
        not_data = []
        newline = {}
        for col, val in zip(colheaders, vals):
            if col in data["data_cols"]:
                try:
                    val = float(val)
                except ValueError:
                    print("value " + val + " of col " + col + " is not data.")
                    not_data += [col]
            newline[col] = val
        if len(not_data) == len(data["data_cols"]):
            print("it looks like there is another data set appended!")
            if multiset:
                print("continuing to next set.")
                numerize(data)
                datasets += [data.copy()]
                colheaders = [val.strip() for val in vals]
                data = get_empty_set(
                    set(colheaders), timestamp=timestamp, data_type="SPEC"
                )
                continue
            else:
                print("returning first set.")
                numerize(data)
                return data
        else:
            for col in not_data:
                data["data_cols"].remove(col)
                print("column " + col + " removed from 'data_cols '.")

        for col, val in zip(colheaders, vals):
            data[col] += [newline[col]]

    numerize(data)
    datasets += [data]
    if verbose:
        print("function 'load_from_csv' finished!")
    if multiset:
        return datasets
    return data


def read_macro(file):
    with open(file) as macro:
        lines = macro.readlines()
    lines = remove_comments(lines)
    settings = {
        "tth": [],
        "alpha": [],
        "savepath": [],
        "newfile": [],
        "measurements": [],
    }
    for line in lines:
        # print(line)
        tth_match = re.search("umv tth " + float_match, line)
        if tth_match:
            # print('got tth!')
            settings["tth"] += [float(tth_match.group()[7:])]
            continue
        alpha_match = re.search("umv th " + float_match, line)
        if alpha_match:
            settings["alpha"] += [float(alpha_match.group()[6:])]
            continue
        if "pd savepath" in line:
            settings["savepath"] += [line[12:]]
            continue
        if "newfile " in line:
            settings["newfile"] += [line[8:]]
            continue
        if "_timescan " in line or "ascan " in line or "pdascan " in line:
            settings["measurements"] += [line]
            continue
    return settings
