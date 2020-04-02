#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 23:10:02 2020

@author: scott
"""


def fix_timecols(data):
    if "timecols" not in data:
        return

    new_timecols = {}
    for col, tcol in data["timecols"].items():
        if col[-2:] == "-x" and tcol[-2:] == "-y":
            new_timecols[tcol] = col
        else:
            new_timecols[col] = tcol

    data["timecols"] = new_timecols
    return data
