# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 21:36:53 2017

This module will capture all the modelling from Chapter 3.1 in Scott's MSc 
thesis for a range of electrolytes, as well as functions interfacing with
EC_MS data to enable predictions of pH change etc.

@author: soren
"""

import os
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import odeint
from matplotlib import pyplot as plt

from .Object_Files import lines_to_dictionary
from . import Chem

data_directory = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data"
with open(data_directory + os.sep + "Electrolytes.txt") as f:
    electrolyte_lines = f.readlines()
    electrolyte_dict = lines_to_dictionary(electrolyte_lines)

Kw = 1.0e-14


def get_electrolyte_type(cation, anion):
    if anion in ["CO3--", "HCO3-", "H2CO3", "CO2", "CO3", "HCO3"]:
        return "carbonate"
    if anion in ["PO4---", "HPO4--", "H2PO4-", "H3PO4", "PO4", "HPO4", "H2PO4"]:
        return "phosphate"
    if anion in ["SO4--", "HSO4-", "H2SO4", "SO4", "HSO4"]:
        return "sulfate"
    if cation in ["NH4+", "NH4"]:
        return "ammonium"
    if cation in ["H+", "H"]:
        return "strong acid"
    if anion in ["OH-", "OH"]:
        return "strong base"
    if anion is not None and cation is not None:
        print("Electrolyte assumed to be neutral salt.")
        return "salt"
    print("Error: unknown electrolyte type")


def read_charge(ion):
    z = 0
    while ion[-1] == "+":
        z += 1
        ion = ion[:-1]
    while ion[-1] == "-":
        z -= 1
        ion = ion[:-1]
    return z


class Electrolyte:
    def __init__(
        self,
        name=None,
        electrolyte_type=None,
        s=None,
        concentration=None,
        pH=None,
        cation=None,
        anion=None,
        spectator=None,
        E=None,
        F=None,
        p=None,
        verbose=True,
    ):
        """
        Allows for a lot of ways of initializing the electrolyte.
        As in Scott's MSc thesis, E is total acid+anion species, p is
        equilibrium pressure of gas.
        """

        if electrolyte_type is not None:
            pass
        elif name is not None:
            try:
                electrolyte_type = electrolyte_dict[name]["electrolyte_type"]
            except KeyError:  # then the first argument is an electrolyte type, not the
                # name of a pre-specified electrolyte
                electrolyte_type = name
                name = None
        else:
            electrolyte_type = get_electrolyte_type(cation, anion)
            name = electrolyte_type
        constants = electrolyte_dict[electrolyte_type]

        self.name = name
        self.constants = constants
        self.electrolyte_type = electrolyte_type
        self.cation = cation
        self.spectator = spectator
        self.anion = anion
        self.concentration = concentration
        self.s = s
        self.pH = pH
        self.p = p
        self.E = E
        self.F = F
        self.verbose = verbose

        if name is not None:  # if it's a named electrolyte
            specs = electrolyte_dict[name]
            for (key, value) in specs.items():
                if not hasattr(self, key):
                    setattr(self, key, value)
                elif getattr(self, key) is None:
                    setattr(self, key, value)

        if self.s is None and hasattr(self, "concentration"):
            self.s = self.concentration
        elif self.concentration is None:
            self.concentration = self.s

        self.set_buffer()
        self.equilibrate()

        if self.name is None:
            self.name = (
                str(self.concentration)
                + " M "
                + self.electrolyte_type
                + ", pH = "
                + str(self.pH)
            )
        if self.verbose:
            print("Initialized electolyte: " + self.name + "\n")

    def set_buffer(self):

        self.species = ["H+", "OH-"]

        for (key, value) in self.constants.items():
            if not hasattr(self, key):
                setattr(self, key, value)
            elif getattr(self, key) is None:
                setattr(self, key, value)

        if "buffer" in self.constants:

            self.pKa = np.array(self.constants["pKa"])
            self.Ka = np.power(10, -self.pKa)
            self.ions = self.constants[
                "buffer"
            ]  # list of the species participating in the
            # acid-base equilibrium from most to least protonated
            # self.cation = self.constants['cation']
            # self.spectator = self.constants['spectator']
        else:
            self.ions = []
            if self.spectator is None:
                if "acid" in self.electrolyte_type:
                    self.spectator = self.anion
                else:
                    self.spectator = self.cation

        self.N_buf = len(self.ions)
        self.z = np.array([read_charge(ion) for ion in self.ions])
        self.species = ["H+", "OH-"] + self.ions + [self.spectator]
        self.charge = dict([(ion, read_charge(ion)) for ion in self.species])

        if "Keq" in self.constants:
            self.dissolved = self.constants["gas"] + "(aq)"
            self.species += [self.dissolved]
        if "Kh" in self.constants:
            self.gas = self.constants["gas"] + "(g)"
            self.species += [self.gas]

    def equilibrate(self, **kwargs):
        """
        Updates all other parameters from any two specified parameters.
        Was a bit tricky to code, but I like what I came up with.
        """
        if "buffer" not in self.constants:
            return self.equilibrate_simple(**kwargs)

        variables = ["s"]  # to be in order of decreasing likelihood to be held constant
        if "Kh" in self.constants:  # then we've got equilibrium with a vapor
            variables += ["p"]
        if (
            "Keq" in self.constants
        ):  # then we've got equilibrium with a dissolved species
            variables += ["E"]

        variables += ["F", "pH"]
        fundict = {
            ("s", "p", "pH"): self.pH_sp,
            ("s", "E", "pH"): self.pH_sE,
            ("s", "F", "pH"): self.pH_sF,
            ("s", "pH", "F"): self.F_spH,
            ("s", "pH", "E"): self.E_spH,
            ("s", "pH", "p"): self.p_spH,
            ("F", "pH", "s"): self.s_FpH,
        }

        done = []
        if len(kwargs) > 2:
            print(
                "Error, can't equilibrate: not enough degrees of freedom!\n"
                + " got the following: "
                + str(done)
            )
            return
        elif len(kwargs) > 0:  # then use those two parameters to set the rest:
            for (var, value) in kwargs.items():
                setattr(self, var, value)
                done += [var]
        while len(done) < 2:
            for var in variables:
                if getattr(self, var) is not None and var not in done:
                    done += [var]
                    # print(str(done))
                    break
            else:
                print(
                    "Error, can't equilibrate: Too many degrees of freedom!\n"
                    + " got the following: "
                    + str(done)
                )
                return

        while len(done) < len(variables):
            try:
                (names, fun) = [
                    (n, f)
                    for (n, f) in fundict.items()
                    if n[0] in done and n[1] in done and n[2] not in done
                ][0]
            except IndexError:
                print(
                    "Error, can't equilibrate: no function appropriate to use.\n"
                    + " got the following: "
                    + str(done)
                )
                return
            print("calling function " + str(fun))  # debugging
            fun()
            #            print('names = ' + str(names) + ', vals = ' +
            #                str([getattr(self, names[0]), getattr(self, names[1])]) + ', a = ' + str(fun()))
            done += [names[2]]

        if self.verbose:
            print("Equilibrated successfully! pH = " + str(self.pH))

    def equilibrate_simple(self, **kwargs):
        variables = ["s", "pH"]
        fundict = {"s": self.pH_s, "pH": self.s_pH}
        if len(kwargs) > 1:
            print(
                "Error, can't equilibrate: not enough degrees of freedom!\n"
                + " got the following: "
                + str(kwargs.keys())
            )
            return
        elif len(kwargs) == 1:  # then use those two parameters to set the rest:
            for var, value in kwargs.items():
                setattr(self, var, value)
                fun = fundict[var]
        else:
            for var in variables:
                if getattr(self, var) is not None:
                    fun = fundict[var]
                    break
            else:
                print(
                    "Error, can't equilibrate: Too many degrees of freedom!\n"
                    + " simply got nothing."
                )
                return
        fun()

        if self.verbose:
            print("Equilibrated successfully! pH = " + str(self.pH))

    def get_concentrations(self, **kwargs):
        if len(kwargs) > 0:
            self.equilibrate(**kwargs)

        bufvec = self.buffer_vec(self.pH)
        bv = sum(bufvec)
        self.conc = {}
        for b, ion in zip(bufvec, self.ions):
            self.conc[ion] = self.F * b / bv

        self.conc["H+"] = np.power(10.0, -self.pH)
        self.conc["OH-"] = Kw / self.conc["H+"]
        if hasattr(self, "dissolved"):
            self.conc[self.dissolved] = self.E - self.F
        if hasattr(self, "gas"):
            self.conc[self.gas] = self.conc[self.dissolved] * self.constants["Kh"]
        self.conc[self.spectator] = self.s

        return self.conc

    def get_conductivity(self, **kwargs):
        if len(kwargs) > 0 or not hasattr(self, "conc"):
            self.get_concentrations(**kwargs)
        self.kappa = 0
        self.kap = {}
        for ion, conc in self.conc.items():
            #         print(ion)
            charge = read_charge(ion)
            if charge != 0:
                self.kap[ion] = (
                    conc
                    * 1e3
                    * Chem.Far
                    * np.abs(charge)
                    * electrolyte_dict["mobility"][ion]
                )
                # units: mol/l * l/m^3 * C/mol * m^2/(V*s) = A/(V*m) = S/m
                self.kappa += self.kap[ion]
        return self.kappa

    def get_conductivities(self, **kwargs):
        self.get_conductivity(**kwargs)
        return self.kap

    def pHsolve(self, residual):
        """
        Whatever equilibrium model I use, this is how the actual numerical 
        solution will go down. Residual is the total net charge concentration 
        as a function of pH, which of course gives 0 for the right pH.
        """
        if self.pH is not None:
            pH0 = self.pH
        else:
            pH0 = 7
        try:
            pH = brentq(residual, pH0 - 0.1, pH0 + 0.1)
        except ValueError:
            if self.verbose:
                print("jump from pH = " + str(pH0) + " : brentq from scratch.")
            pH = brentq(residual, -6, 20)
        self.pH = pH
        return pH

    def pH_s(self, s=None):
        """
        sets self.pH using s-F equilibrium, i.e. excluding dissolved gas
        """
        if s is None:
            s = self.s

        sc = s * self.charge[self.spectator]  # spectator charge

        def residual(pH):
            return sc + np.power(10, -pH) - Kw * np.power(10, pH)

        return self.pHsolve(residual)

    def s_pH(self, pH=None):

        if pH is None:
            pH = self.pH

        x = np.power(10.0, -pH)

        sc = Kw / x - x  # spectator charge
        s = sc / self.charge[self.spectator]
        self.s = s
        return s

    def buffer_vec(self, pH):
        """
        Useful concept for acid-base equilibrium with multiple Ka's:
        bufvec = [1, Ka1/x, Ka2*Ka1/x^2, ...]
        where x is proton concentration, and pKa1<pKa2<...
        """
        return np.array(
            [
                np.product(self.Ka[0:i]) / np.power(10.0, -pH * i)
                for i in range(self.N_buf)
            ]
        )
        # if 10 rather than 10.0, I get nans

    def pH_sp(self, s=None, p=None):
        """
        sets self.pH using s-p equilibrium. 
        Only works for electrolytes in eq. with gas, e.g. carbonate
        """
        if s is None:
            s = self.s
        if p is None:
            p = self.p
        Keq = self.constants["Keq"]
        Kh = self.constants["Kh"]

        def buffer_charge(pH):
            bufvec = self.buffer_vec(pH)
            return p * Keq / Kh * np.dot(self.z, bufvec)

        sc = s * self.charge[self.spectator]  # spectator charge

        def residual(pH):
            return sc + np.power(10, -pH) - Kw * np.power(10, pH) + buffer_charge(pH)

        return self.pHsolve(residual)

    def pH_sE(self, s=None, E=None):
        """
        sets self.pH using s-E equilibrium, i.e. including dissolved gas
        """
        if s is None:
            s = self.s
        if E is None:
            E = self.E
        try:
            Keq = self.constants["Keq"]
        except KeyError:
            if self.verbose:
                print("no gas equilibrium in electrolyte " + self.name)
            return self.pH_sF(s, E)

        def buffer_charge(pH):
            bufvec = self.buffer_vec(pH)
            return E * np.dot(self.z, bufvec) / (1 / Keq + np.sum(bufvec))

        sc = s * self.charge[self.spectator]  # spectator charge

        def residual(pH):
            return sc + np.power(10, -pH) - Kw * np.power(10, pH) + buffer_charge(pH)

        return self.pHsolve(residual)

    def pH_sF(self, s=None, F=None):
        """
        sets self.pH using s-F equilibrium, i.e. excluding dissolved gas
        """
        if s is None:
            s = self.s
        if F is None:
            F = self.F

        def buffer_charge(pH):
            bufvec = self.buffer_vec(pH)
            return F * np.dot(self.z, bufvec) / np.sum(bufvec)

        sc = s * self.charge[self.spectator]  # spectator charge

        def residual(pH):
            return sc + np.power(10, -pH) - Kw * np.power(10, pH) + buffer_charge(pH)

        return self.pHsolve(residual)

    def p_spH(self, s=None, pH=None):
        if s is None:
            s = self.s
        if pH is None:
            pH = self.pH
        try:
            Kh = self.constants["Kh"]
            Keq = self.constants["Keq"]
        except KeyError:
            # print('no gas equilibrium in electrolyte ' + self.name)
            return self.E_spH(s, pH)
        #        print('pH = ' + str(pH) + ', s = ' + str(s))
        x = np.power(10.0, -pH)  # 10 instead of 10.0 gives error. wtf numpy?
        bufvec = self.buffer_vec(pH)
        #       print('bufvec = ' + str(bufvec) + '\ndot product = ' + str(np.dot(bufvec, self.z)))
        sc = s * self.charge[self.spectator]  # spectator charge
        p = (-sc - x + Kw / x) * Kh / (Keq * np.dot(bufvec, self.z))
        self.p = p
        return p

    def E_spH(self, s=None, pH=None):
        if s is None:
            s = self.s
        if pH is None:
            pH = self.pH
        try:
            Keq = self.constants["Keq"]
        except KeyError:
            # print('no gas equilibrium in electrolyte ' + self.name)
            return self.F_spH(s, pH)
        x = np.power(10.0, -pH)
        bufvec = self.buffer_vec(pH)
        sc = s * self.charge[self.spectator]  # spectator charge
        E = (-sc - x + Kw / x) * (sum(bufvec) + 1 / Keq) / np.dot(bufvec, self.z)
        self.E = E
        return E

    def F_spH(self, s=None, pH=None):
        if s is None:
            s = self.s
        if pH is None:
            pH = self.pH
        x = np.power(10.0, -pH)
        bufvec = self.buffer_vec(pH)
        sc = s * self.charge[self.spectator]  # spectator charge
        F = (-sc - x + Kw / x) * sum(bufvec) / np.dot(bufvec, self.z)
        self.F = F
        return F

    def s_FpH(self, F=None, pH=None):
        if F is None:
            F = self.F
        if pH is None:
            pH = self.pH

        x = np.power(10.0, -pH)
        bufvec = self.buffer_vec(pH)
        sc = Kw / x - x - F * np.dot(bufvec, self.z) / sum(bufvec)  # spectator charge
        s = sc / self.charge[self.spectator]
        self.s = s
        return s


def electrolysis_ode(quantities, t, pars):
    """
    It'll take some work to make this completely general. Designed now for s-E equilibrium
    for carbonate.
    Also this is in serious need of vectorization.
    """
    electrolyte = pars[0]
    equilibrium_type = pars[1]  # for example: sE, or sF
    cq_dot = pars[2](t)

    quantitiesdict = {
        equilibrium_type[0]: quantities[0],
        equilibrium_type[1]: quantities[1],
    }
    # quantitiesdict is a **kwarg dictionary for inputting the parameters specified
    # by the equilibrium type.
    kap = electrolyte.get_conductivities(**quantitiesdict)
    # this will equilibrate the electrolyte. Then I can just read the conductivities
    kappa = electrolyte.kappa
    ddt = {}
    ddt["F"] = 0
    for ion, k in kap.items():
        #        print('ion: ' + ion)
        if ion == electrolyte.spectator:
            ddt["s"] = k / (read_charge(ion) * kappa) * cq_dot
        elif ion in electrolyte.ions:
            ddt["F"] += k / (read_charge(ion) * kappa) * cq_dot
    ddt["E"] = ddt["F"]

    dquantitiesdt = np.array([ddt[q] for q in list(equilibrium_type)])

    return dquantitiesdt


def electrolysis_simple_ode(s, t, pars):
    electrolyte = pars[0]
    cq_dot = pars[1](t)

    # quantitiesdict is a **kwarg dictionary for inputting the parameters specified
    # by the equilibrium type.
    kap = electrolyte.get_conductivities(s=s)
    k = kap[electrolyte.spectator]
    z = electrolyte.charge[electrolyte.spectator]
    # this will equilibrate the electrolyte. Then I can just read the conductivities
    kappa = electrolyte.kappa

    dsdt = k / (z * kappa) * cq_dot

    return dsdt


def electrolyze(
    electrolyte="standard",
    equilibrium_type="sE",
    tj=None,
    tpulse=60,
    tspan=None,
    colors={},
    j_el=-5,
    L=100e-6,
    ax="new",
    leg=True,
):
    """
    calculates pH and acid/base concentrations as a function of time during 
    electrolysis, starting with a given electrolyte state and assuming
    proton-consuming (j_el<0) or proton-releasing (j_el>0) reactions and an 
    isolated film (diffusion layer) of thickness L over the electrode.
    """
    if tj is None:
        if tspan is None:
            tspan = [0, tpulse]
        tvec = np.linspace(tspan[0], tspan[-1], 1000)

        def j_fun(t):
            if t < 0:
                return 0
            if t < tpulse:
                return j_el
            return 0

    else:
        t_in = tj[0]
        j_in = tj[1]
        tvec = t_in

        def j_fun(t):
            if t < t_in[0]:
                return 0
            if t < t_in[-1]:
                return np.interp(t, t_in, j_in)
            return 0

    def cq_dot_fun(t):
        return -j_fun(t) / (Chem.Far * L) * 1e-2
        # units: mA/cm^2 / (C/mol * m) * (A/m^2/(mA/cm^2)) (M / (mol/m^3)) = M / s

    if type(electrolyte) is str:
        electrolyte = Electrolyte(electrolyte)

    electrolyte.equilibrate()
    electrolyte.verbose = False  # otherwise it really screams.

    vec = {}
    quantitiesdict = []
    if "buffer" in electrolyte.constants:
        pars = (electrolyte, equilibrium_type, cq_dot_fun)
        var_0 = equilibrium_type[0]
        var_1 = equilibrium_type[1]

        quantities0 = [getattr(electrolyte, var_0), getattr(electrolyte, var_1)]
        solution = odeint(electrolysis_ode, quantities0, tvec, args=(pars,))
        vec[var_0] = solution[:, 0]
        vec[var_1] = solution[:, 1]
        for (q0, q1) in zip(vec[var_0], vec[var_1]):
            quantitiesdict += [{var_0: q0, var_1: q1}]
    else:
        s0 = electrolyte.s
        pars = (electrolyte, cq_dot_fun)
        solution = odeint(electrolysis_simple_ode, s0, tvec, args=(pars,))
        vec["s"] = solution
        quantitiesdict = [{"s": q} for q in vec["s"]]

    pHvec = []
    y = {}
    for species in electrolyte.get_concentrations().keys():
        if "(g)" not in species:
            print("preparing space to get: " + species)
            y[species] = []
    for qdict in quantitiesdict:
        for species, conc in electrolyte.get_concentrations(**qdict).items():
            #    print('geting data for: ' + species)
            if species in y:
                y[species] += [conc]
        pHvec += [electrolyte.pH]

    default_colors = ["m", "b", "c", "g", "y", "r", "0.5"]

    if ax == "new":
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax2 = ax1.twinx()
    for species, vec in y.items():
        if species in colors:
            color = colors[species]
        else:
            color = default_colors.pop()
        print("plotting: " + species)
        ax1.plot(tvec, vec, color=color, label=species)
    ax2.plot(tvec, pHvec, color="k", label="pH")
    if leg:
        ax1.legend()
    ax1.set_ylabel("concentration / M")
    ax1.set_xlabel("time / s")
    ax2.set_ylabel("pH")

    return [ax1, ax2]


if __name__ == "__main__":
    plt.close("all")
    el1 = Electrolyte("standard", verbose=True)
    print("saturating electrolyte 1 with 1 bar gas...")
    el1.equilibrate(p=1)
    print("saturated, pH = " + str(el1.pH))

    el2 = Electrolyte(
        cation="K", anion="PO4---", concentration=1.0, pH=12, verbose=False
    )
    titrateit = False
    if titrateit:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax2 = ax1.twinx()
        F = 0.5
        el2.F = F
        y = {}
        pHvec = []
        for species in el2.get_concentrations().keys():
            if species in el2.ions:
                y[species] = []
        svec = np.linspace(0, 2, 1000)
        for s in svec:
            for species, conc in el2.get_concentrations(s=s, F=F).items():
                if species in y:
                    y[species] += [conc]
            pHvec += [el2.pH]

        colors = ["m", "b", "c", "g", "y", "r"]
        for (species, vec), color in zip(y.items(), colors):
            ax1.plot(svec, vec, color=color, label=species)
        ax2.plot(svec, pHvec, color="k", label="pH")
        ax1.legend()
        ax1.set_ylabel("concentration / M")
        ax1.set_xlabel("K+ concentration / M")
        ax2.set_ylabel("pH")
        ax1.set_title("total acid/base concentration = " + str(F) + " M")

    electrolyzeit = True
    if electrolyzeit:
        el3 = Electrolyte("phosphate", F=0.5, s=0.05)
        electrolyze(el3, equilibrium_type="sF", tpulse=100, j_el=-10, ax="new")

    el4 = Electrolyte("perchloric acid", s=0.1)

    el6 = Electrolyte("hydroxide", s=1)

    electrolyze(el6, tpulse=100, j_el=10, ax="new")

    # print(el1.pH_sp(s=0.1, p=400e-6))
