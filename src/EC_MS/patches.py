#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 23:10:02 2020

@author: scott
"""


def fix_timecols(data):
    """Update and fix old errors in the timecols item of old data dictionaries"""

    if "timecols" not in data:
        return

    new_timecols = {}

    # some old pickles have timecols as tuples:
    if not isinstance(data["timecols"], dict):
        data["timecols"] = dict(data["timecols"])

    for col, tcol in data["timecols"].items():
        if col.endswith("-x") and tcol.endswith("-y"):
            new_timecols[tcol] = col
        else:
            new_timecols[col] = tcol

    data["timecols"] = new_timecols
    return data
