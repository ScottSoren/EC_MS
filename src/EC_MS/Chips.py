#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:26:33 2017

@author: scott
"""


from __future__ import print_function, division
import os
import numpy as np
from matplotlib import pyplot as plt

from . import Chem
from .Molecules import Molecule
from .Object_Files import structure_to_lines, lines_to_dictionary
from .Object_Files import lines_to_structure, date_scott, update_lines

data_directory = (
    os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "chips"
)
cwd = os.getcwd()
# for python2:
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


design_parameters = {}


design_parameters["membrane_chip_02"] = dict(
    w_cap=6e-6,  # capillary width / [m]
    h_cap=6e-6,  # capillary height / [m]
    l_cap=1e-3,  # capillary length / [m]
    p=1e5,  # chip pressure /
    T=298,
)


class Chip:
    def __init__(
        self, name="unknown", chip_type="membrane_chip_02", updatefile=False, **kwargs
    ):
        """
        Chip parameters taken from Daniel's thesis
        """
        self.parameters = {}

        self.name = name

        if not self.reset():  # reset returns false if there is no existing data file
            try:
                self.parameters = design_parameters[chip_type]
            except KeyError:
                print("WARNING!!! no parameters available for chip '" + name + "'.")
            else:
                self.write(name=name)

        parameters = self.parameters
        parameters.update(kwargs)
        for key, value in parameters.items():
            setattr(self, key, value)
            if updatefile:
                self.write((key, value))

    def capillary_flow(
        self, gas="He", w_cap=None, h_cap=None, l_cap=None, T=None, p=None
    ):
        """
        adapted from Membrane_chip.py,
        equations from Daniel Trimarco's PhD Thesis

        Returns the flux in molecules/s of a carrier gas through the
        chip capillary.
        As the flow starts out viscous, at low analyte production rates,
        an analyte flux is simply the analyte's mol fraction in the chip times
        this flux.

        We assume that flow is governed by the bulk properties (viscosity)
        of the carrier gas and the molecular properties (diameter, mass)
        of the analyte.
        """

        if type(gas) is str:
            gas = Molecule(gas)

        """ #This does not work! I do not know why.
        print(self.w_cap)
        for var in 'w_cap', 'h_cap', 'l_cap', 'T', 'p':
            if locals()[var] is None:
                locals()[var] = getattr(self, var)
        """

        # This is actually the most compact way I can figure out to do this.
        w = next(x for x in [w_cap, self.w_cap] if x is not None)
        h = next(x for x in [h_cap, self.h_cap] if x is not None)
        l = next(x for x in [l_cap, self.l_cap] if x is not None)
        T = next(x for x in [T, self.T] if x is not None)
        p = next(x for x in [p, self.p] if x is not None)

        s = gas.molecule_diameter  # molecular diameter in m
        m = gas.molecule_mass  # molecular mass in kg
        eta = gas.dynamic_viscosity  # viscosity in Pa*s

        d = ((w * h) / np.pi) ** 0.5 * 2  # hydraulic diameter
        # d=4.4e-6  #used in Henriksen2009
        a = d / 2  # hydraulic radius
        p_1 = p  # pressure at start of capillary (chip pressure)
        lam = d  # mean free path of transition from visc. to mol. flow

        p_t = Chem.kB * T / (2 ** 0.5 * np.pi * s ** 2 * lam)  # transition pressure

        p_2 = 0  # pressure at end of capillary (vacuum)

        p_m = (p_1 + p_t) / 2  # average pressure in the visc. + trans. flow region

        v_m = (8 * Chem.kB * T / (np.pi * m)) ** 0.5  # mean molecular velocity

        nu = (m / (Chem.kB * T)) ** 0.5  # a resiprocal velocity used for short-hand
        # dumb, but inherited straight from Henrikson2009 through all of
        # Daniel's work up to his thesis, where it was finally dropped, which
        # unfortunatly makes that term of the equation a bit silly-looking.

        N_dot = (
            1
            / (Chem.kB * T)
            * 1
            / l
            * (
                (
                    a ** 4 * np.pi / (8 * eta) * p_m
                    + a ** 3
                    * 2
                    * np.pi
                    / 3
                    * v_m
                    * (1 + 2 * a * nu * p_m / eta)
                    / (1 + 2.48 * a * nu * p_m / eta)
                )
                * (p_1 - p_t)
                + a ** 3 * 2 * np.pi / 3 * v_m * (p_t - p_2)
            )
        )
        return N_dot

    def write(self, a=None, attr=None, name=None, *args, **kwargs):
        """
        Identical to the function in class Molecule of same name
        ... could move this to Object_Files, but a bit tricky

        this is supposed to be a versitile tool for writing to the Molecule's
        data file. Whether the added intricacy will be worth the lines of code
        it saves, I can't tell yet.

        Writes in one of the following ways:

        1. If the name of an attribute is given, that attribute is written.
        2. If a is a string, simply write a to the molecule's datafile.
        3. If a is a function, then it's a function that writes what you want
           to the data file given in the first argument
        4. If a is a dictionary, list, or tuple, convert it to lines according to
           the system encoded in Object_Files.
        5. If a is not given but keyword arguments are, write **kwargs to the
           data file according to the system encoded in Object_Files.
        """

        if name is None:
            name = self.name
        if attr is not None:
            a = (attr, getattr(self, attr))
        elif a is None:
            if len(kwargs) == 0:
                print("nothing to write.")
                return
            else:
                a = kwargs

        cwd = os.getcwd()
        os.chdir(data_directory)
        file_name = name + ".txt"

        with open(file_name, "a") as f:
            if callable(a):
                a(f, *args, **kwargs)
            elif type(a) is str:
                if a[-1] != "\n":
                    a += "\n"
                f.write(a)
            elif type(a) in (list, dict):
                if "key" in kwargs.keys():
                    lines = structure_to_lines(a, preamble=kwargs["key"])
                else:
                    lines = structure_to_lines(a)
                f.writelines(lines)
            elif type(a) is tuple:
                lines = structure_to_lines(a[1], preamble=a[0])
                f.writelines(lines)
            else:
                print("Couldn" "t write " + str(a))
        os.chdir(cwd)

    def save(self, name):
        """
        Save the important stuff in the chip, namely
        """
        if name is None:
            name = self.name
        else:
            self.name = name
        print("Saving info on Chip " + name + " to file")
        file_name = name + ".txt"
        cwd = os.getcwd()
        os.chdir(data_directory)
        with open(file_name, "w") as f:
            f.write("# created " + date_scott() + "\n")
        os.chdir(cwd)
        for key, value in self.parameters.items():
            setattr(self, key, value)
            self.write((key, value))

    def reset(self):
        file_name = self.name + ".txt"
        try:
            cwd = os.getcwd()
            os.chdir(data_directory)
            f = open(file_name, "r")
        except FileNotFoundError:
            print("no file for " + self.name)
            os.chdir(cwd)
            return False

        print("loaded info on Chip " + self.name + " from file.")
        lines = f.readlines()
        f.close()
        os.chdir(cwd)

        dictionary = lines_to_dictionary(lines)
        self.parameters.update(dictionary)
        for key, value in dictionary.items():
            setattr(self, key, value)
        return True


if __name__ == "__main__":

    c = Chip(name="SI-3iv1-1-C5")

    N_dot_He = c.capillary_flow("He")
    print("He flow is " + str(N_dot_He / Chem.NA * 1e9) + " nmol/s")
    N_dot_CO = c.capillary_flow("CO")
    print("CO flow is " + str(N_dot_CO / Chem.NA * 1e9) + " nmol/s")
    N_dot_CH4 = c.capillary_flow("CH4")
    print("CH4 flow is " + str(N_dot_CH4 / Chem.NA * 1e9) + " nmol/s")
