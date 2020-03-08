# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 11:06:17 2020

@author: Spectro Inlets Lab_1
"""


import numpy as np

from .Combining import remove_negatives


def read_PVMassSpec(file_path, delim="\t"):
    with open(file_path, "r") as f:
        lines = f.readlines()

    dataset = {"file": file_path, "header": "", "timecols": {}}
    N_col_head = len(
        lines
    )  # this will decrease when the loop knows when the column header line is comming

    nondata_cols = []  # this will store abstime, so I don't have to parse.
    for n, line in enumerate(lines):
        l = line.strip()
        if n < N_col_head:
            dataset["header"] = dataset["header"].join(line)
            if len(l) == 0:
                N_col_head = n + 1
        elif n == N_col_head:
            data_cols = l.split(delim)
            for col in data_cols:
                dataset[col] = np.array([])
                dataset["timecols"][col] = "Time Relative (sec)"
        elif n > N_col_head:
            for col, val in zip(data_cols, l.split(delim)):
                if col in nondata_cols:
                    dataset[col] += [val]
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
                    dataset[col] = list(dataset[col])
                    dataset[col] += [val]
                    nondata_cols += [col]
                else:
                    dataset[col] = np.append(dataset[col], x)

    dataset["data_cols"] = set([col for col in data_cols if col not in nondata_cols])

    return dataset


def rename_PVMassSpec_cols(dataset):
    data_cols = dataset["data_cols"].copy()
    timecol = "Time Relative (sec)"
    for col in data_cols:
        if "_amu" in col:
            mass = "M" + col.split("_amu")[0]
            xcol = mass + "-x"
            ycol = mass + "-y"
            dataset[ycol] = dataset[col]
            dataset[xcol] = dataset[timecol]
            dataset["data_cols"].add(xcol)
            dataset["data_cols"].add(ycol)
            dataset["timecols"][ycol] = xcol
            dataset["timecols"][xcol] = xcol
    # remove_negatives(dataset) # not sure if this is useful.
    return dataset  # not really necessary
