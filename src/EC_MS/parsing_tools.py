# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:45:37 2020

@author: scott
"""

import re
import time, datetime
import numpy as np
import os

"""
"""


float_match = "[-]?\d+[\.]?\d*(e[-]?\d+)?"  # matches floats like '-3.5e4' or '7' or '245.13' or '1e-15'
# note, no white space included on the ends! Seems to work fine.
timestamp_match = "([0-9]{2}[:_.]){2}[0-9]{2}"  # matches timestamps like '14:23:01' or '15_42_15' or '09.56.25'
date_match = "([0-9]{2}[/\-\.]){2}[0-9]{4}"  # matches dates like '01/15/2018' or '09-07-2016' or '04.20.2019'
date_match_2 = (
    "[0-9]{4}([/-][0-9]{2}){2}"  # matches dates like '2018/01/15' or '2018-09-07'
)


def numerize(data):
    for col in data["data_cols"]:  # numerize!
        data[col] = np.array(data[col])


def get_empty_set(cols, **kwargs):
    # get the colheaders and make space for data
    # cols should be a set
    data = {}
    data.update(kwargs)
    for col in cols:
        data[col] = []
    data["data_cols"] = cols
    return data


def parse_timezone(tz=None):
    """
    Gets a timezone object from a timezone string. Includes some abbreviations
    useful for Scott. If the input is not a string, it is returned as is.
    """
    abbreviations = {
        "CA": "US/Pacific",
        "DK": "Europe/Copenhagen",
    }
    if tz is not None:
        import pytz  # all three seem necessary for dealing with timezones.
    if tz in abbreviations:
        return pytz.timezone(abbreviations[tz])
    elif type(tz) is str:
        return pytz.timezone(tz)
    return tz


def parse_date(line):
    # ^  mm/dd/yyyy, as EC lab does
    # older EC-Lab seems to have dashes in date, and newer has slashes.
    # Both seem to save month before day, regardless of where the data was taken or .mpt exported.
    d1 = re.search(date_match, line)
    if d1:
        date1 = d1.group()
        yyyy, mm, dd = date1[-4:], date1[:2], date1[3:5]
        date = yyyy + "/" + mm + "/" + dd
        return date

    # ^  yyyy/mm/dd, as cinfdata does
    d2 = re.search(date_match_2, line)
    if d2:
        date2 = d2.group()
        yyyy, mm, dd = date2[:4], date2[5:7], date2[-2:]
        date = yyyy + "/" + mm + "/" + dd
        return date

    # shit, some dates are written in long form. This'll be a tough RE exercise. Here we go!
    # I want to match dates like 'Apr. 12, 2019' or 'August 1, 2019'
    month_names = {
        "01": "January",
        "02": "February",
        "03": "March",
        "04": "April",
        "05": "May",
        "06": "June",
        "07": "July",
        "08": "August",
        "09": "September",
        "10": "October",
        "11": "November",
        "12": "December",
    }
    month_match = (
        "("
        + "".join(list([v[:3] + "|" for v in month_names.values()]))[:-1]
        + ")"
        + "[a-z]*(\.)?"
    )
    date_match_3 = month_match + " [0-9]+, [0-9]{4}"
    d3 = re.search(date_match_3, line)
    if d3:
        date3 = d3.group()
        month = re.search(month_match, date3).group()
        mm = next(key for key, value in month_names.items() if value[:3] == month[:3])
        yyyy = date3[-4:]
        dd = re.search("[0-9]+,", date3).group()[:-1]
        date = yyyy + "/" + mm + "/" + dd
        return date

    print("can't find date in line '" + line + "'. parse_date() is returning None.")
    return None


def timestring_to_epoch_time(
    timestring, date=None, tz=None, verbose=True, form="%Y/%m/%d %H:%M:%S", out="tstamp"
):
    """
    A way to convert a number of timestrings read from my data into a standard-formatted
    date and time, and then to an epoch (unix) time.

    tz is the Timezone, which is only strictly necessary when synchronizing
    data taken at different places or accross dst at a place with different dst
    implementation than the local (dst is stupid!).
    If timezone is a number, it is interpreted as the offset from GMT of the
    data.

    The epoch time is referred to here and elsewhere as tstamp.
    """
    if tz is not None:
        import pytz  # all three seem necessary for dealing with timezones.

        tz = parse_timezone(tz)
        if verbose:
            print("getting epoch time given a timestamp local to " + str(tz))
        epoch = pytz.utc.localize(datetime.datetime.utcfromtimestamp(0))

    if timestring == "now":
        return time.time()
    elif type(timestring) is time.struct_time:
        if verbose:
            print(
                "'timestring_to_epoch_time' revieved a time.struct_time object. "
                + "Returning the corresponding epoch time."
            )
        return time.mktime(timestring)

    elif type(timestring) is not str:
        if verbose:
            print(
                "WARNING: 'timestamp_to_unix_time' didn't receive a string. "
                + "Received: "
                + str(timestring)
                + " . Returning the argument."
            )
        return timestring

    if len(timestring) > 8:
        try:
            # print(timestring) # debugging
            timestamp = re.search(timestamp_match, timestring).group()
            hh = int(timestamp[0:2])
            if "PM" in timestring and not hh == 12:
                # Holy fuck the whole AM/PM thing is stupid
                timestamp = str(hh + 12) + timestamp[2:]
            elif "AM" in timestring and hh == 12:
                timestamp = "00" + timestamp[2:]
        except AttributeError:
            if verbose:
                print(
                    "WARNING: I got no clue what you're talking 'bout "
                    + "when you say "
                    + timestamp
                    + ". It didn't match \."
                    + timestamp_match
                    + "'. Assuming you want 00:00:00"
                )
            timestamp = "00:00:00"
        if verbose:
            print("found timestamp = " + timestamp)
    else:
        timestamp = timestring
    if "_" in timestamp:
        timestamp = timestamp.replace(
            "_", ":"
        )  # otherwise time.strptime below crashes.
    if "." in timestamp:
        timestamp = timestamp.replace(
            ".", ":"
        )  # otherwise time.strptime below crashes.

    if date is None:
        if verbose:
            print(
                "'timestring_to_epoch_time' is assuming"
                + " the date is in the timestring."
            )
        date = parse_date(timestring)

    if date is None:
        if verbose:
            print(
                "couldn't figure out the date for " + timestring + ". Assuming today."
            )
        date = "today"
    if date == "today":
        date = time.strftime("%Y/%m/%d")
    # print('timestring = ' + timestring) # debugging
    if tz is None:
        if "-" in date:
            date = date.replace("-", "/")  # 18D08
        struct = time.strptime(date + " " + timestamp, form)
        tstamp = time.mktime(struct)
    else:
        dt_naive = datetime.datetime.strptime(date + " " + timestamp, form)
        dt = tz.localize(dt_naive)
        tstamp = (dt - epoch).total_seconds()

    if out == "all":
        return tstamp, date, timestamp
    return tstamp


def epoch_time_to_timestamp(tstamp, tz=None, verbose=True):
    """
    tz is the Timezone, which is only strictly necessary when synchronizing
    data taken at different places or accross dst at a place with different dst
    implementation than the local (dst is stupid!).
    If timezone is a number, it is interpreted as the offset from GMT of the
    data in hours (+1 for Denmark, -8 for California)
    """
    if tz is None:
        struct = time.localtime(tstamp)
    else:
        tz = parse_timezone(tz)
        if verbose:
            print("getting the timestamp local to " + str(tz) + " from epoch time.")
        dt_utc = datetime.datetime.utcfromtimestamp(tstamp)
        dt_tz = tz.fromutc(dt_utc)
        struct = dt_tz.timetuple()
    hh = str(struct.tm_hour)
    if len(hh) == 1:
        hh = "0" + hh
    mm = str(struct.tm_min)
    if len(hh) == 1:
        mm = "0" + mm
    ss = str(struct.tm_sec)
    if len(hh) == 1:
        ss = "0" + ss
    timestamp = hh + ":" + mm + ":" + ss
    return timestamp


def timetag_to_timestamp(filename):
    """
    Converts a time tag of format _<hh>h<mm>m<ss>_ to timestamp of format
    <hh>:<mm>:<ss>, which is what synchronize reads (I should maybe change
    synchronize to work with the unix epoch timestamp instead...)
    The time tag is something we write in the file names to give an approximate
    sense of when a measurement is started. It can be used as the measurement
    start time if there really is no better way.
    I can't beleive SPEC doesn't save time. I will pressure them to fix this.
    """
    hm_match = re.search(r"_[0-9]{2}h[0-9]{2}", filename)
    hm_str = hm_match.group()
    hh = hm_str[1:3]
    mm = hm_str[-2:]
    ss_match = re.search(r"_[0-9]{2}h[0-9]{2}m[0-9]{2}", filename)
    if ss_match is None:
        ss = "00"
    else:
        ss = ss_match.group()[-2:]
    return hh + ":" + mm + ":" + ss


def get_creation_timestamp(filepath):
    """
    Returns creation timestamp of a file in the format that
    combining.syncrhonize reads.
    The timestamp is local time, not absolute time.
    We need to move to epoch time everywhere!!!
    """
    t = get_creation_time(filepath)
    struct = time.localtime(t)
    hh = str(struct.tm_hour)
    mm = str(struct.tm_minute)
    ss = str(struct.tm_second)
    return hh + ":" + mm + ":" + ss


def get_creation_time(filepath, verbose=True):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    import platform

    if platform.system() == "Windows":
        tstamp = os.path.getctime(filepath)
        if verbose:
            print("In Windows. Using os.path.getctime('" + filepath + "') as tstamp.")
    else:
        stat = os.stat(filepath)
        try:
            tstamp = stat.st_birthtime
            if verbose:
                print(
                    "In linux. Using os.stat('"
                    + filepath
                    + "').st_birthtime as tstamp."
                )
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            tstamp = stat.st_mtime
            if verbose:
                print(
                    "Couldn't get creation time! Returing modified time.\n"
                    + "In linux. Using os.stat('"
                    + filepath
                    + "').st_mtime as tstamp."
                )
    return tstamp


def timestamp_from_file(filepath, verbose=True):
    a = re.search("[0-9]{2}h[0-9]{2}", filepath)
    if a is None:
        if verbose:
            print("trying to read creation time")
        timestamp = get_creation_timestamp(filepath)
    else:
        if verbose:
            print("getting timestamp from filename " + filepath)
        timestamp = timetag_to_timestamp(filepath)
    return timestamp


def remove_comments(lines):
    new_lines = []
    for line in lines:
        if "#" in line:
            line = re.search("^.*\#", line).group()[:-1]
            if re.search(r"\w", line):  # to drop lines that only have comments
                new_lines += [line]
        else:
            new_lines += [line]  # I don't want to get rid of empty lines here
    return new_lines
