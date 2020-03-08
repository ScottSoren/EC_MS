#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 00:21:47 2017

@author: scott


17G30:
There may be something time-saving here:
/home/scott/Dropbox/other_DTU/MSc/Scripts/PYTHON3/Masters Project/15L29 Report Intro
... but I can probably do better from scratch

or, better, build off someone elses work.
how about this?
https://pypi.python.org/pypi/CoolProp/2.2.3
#or a function that looks up from hbcp or something...
googling around does not make it seem easy
"""

import re
import numpy as np
from math import gcd

from .PhysCon import R, Far
from .MolarMasses import get_elements


dfH0 = {  # standard enthalpies of formation / [kJ/mol]
    "H2O(g)": -241.82,
    "H2O(l)": -285.8,
    "CH3CH2OH(g)": -234.8,  # Langes Handbook
    "CH3CH2OH(l)": -277,
    "CH3CH2CH2OH(g)": -256,  # NIST
    "CH3CH2CHO(g)": -188.7,  # NIST
    "C2H6(g)": -84.0,  # NIST
    "HCOOH(g)": -378.6,  # NIST
    "CH4(g)": -74.9,  # NIST
    "CH3CHO(g)": -166.1,
    "CH3OH(g)": -201.2,  # Langes Handbook, John A Dean, 15th edition
    "CO(g)": -110.5,
    "C2H4(g)": 52.4,  # NIST
    "C(s)": 0,
    "O2(g)": 0,
    "H2(g)": 0,
    "Cu2O(s)": -170.71,  # NIST
    "CuO(s)": -156.06,  # NIST
    "Cu(OH)2(s)": -450.37,  # NIST
    "CO2(g)": -393.5,
}
S0 = {  # standard entropy / [J/(mol*K)]
    "H2O(g)": 188.72,
    "H2O(l)": 69.940,
    "CO2(g)": 213.8,
    "CH3CH2OH(g)": 281.6,  # Langes Handbook
    "CH3CH2OH(l)": 161,
    "CH3CH2CH2OH(g)": 322.49,  # NIST
    "CH3CH2CHO(g)": 304.4,  # NIST
    "HCOOH(g)": 248.7,  # NIST
    "CO(g)": 197.7,
    "C2H4(g)": 219.3,  # NIST
    "H2(g)": 130.68,  # NIST
    "C(s)": 5.6,  # NIST
    "N2(g)": 191.61,  # NIST
    "O2(g)": 205.15,  # NIST
    #'C2H6(g)':229.5, #http://bilbo.chm.uri.edu/CHM112/tables/thermtable.htm
    "CH4(g)": 186.7,  # NIST
    "CH3CHO(g)": 263.8,
    "CH3OH(g)": 126.8,  # Langes Handbook, John A Dean, 15th edition
    "Cu(s)": 33.17,  # NIST
    "Cu2O(s)": 92.37,  # NIST
    "CuO(s)": 42.59,  # NIST
    "Cu(OH)2(s)": 108.43,  # NIST
}
dfG0 = {  # standard free energies of formation / [kJ/mol]
    # I need to find a way to querry a reliable database, i.e., NIST.
    # Most of these standard energies are from
    """ Commented out because the source is no longer accessible.
        # bilbo.chm.uri.edu/CHM112/tables/thermtable.htm
        'H2O':-237.1,
        'H2(g)':0,'O2(g)':0, 'Cu':0, 'C(s)':0,
        'H2O(l)':-237.2,'CO2(g)':-394.4,'HCOOH(l)':-346,
        'CO(g)':-137.2,'CH3OH(g)': -162.3,#'CH4(g)':-50.75,
	     'CH3COOH(l)':-389.9,'CH3CHO(g)':-133.4,'CH3CH2OH(g)':-167.9,
        'C2H4(g)':68.12,'C3H8O(g)':-163.0,'C2H2O4(g)':-662.7,
        'CO2(aq)':-386.2,
        'HCOOH(aq)':-356,
        'CH3COOH(aq)':-396.6,
        'C2H2O4(aq)':-697.0, 'CH3OH(l)':-166.4, 'CH3CH2OH(l)':-174.9,
        """
    # http://www2.ucdsb.on.ca/tiss/stretton/database/organic_thermo.htm:
    "C3H6": 62.0,  # 74.7 , ... can't find a reliable one.
    "C3H8": -23.4,
    "C10H22": 33.32,  # https://www.chemeo.com/cid/44-644-8/Decane
    "Cu(s)": 0,
    "Cu2(OH)2CO3(s)": -894.00,  # Kiseleva1992, https://link.springer.com/content/pdf/10.1007%2FBF00204009.pdf
    "e-": 0,
    "H+": 0,  # defining standard state in electrochemistry
    "": 0,  # anticipating the unanticipated
}

# standard pure substance states for formation energies
pure_states = {
    "H": "H2(g)",
    "C": "C(s)",
    "N": "N2(g)",
    "O": "O2(g)",
    "Cu": "Cu(s)",
}

standard_states = dict(
    [(mol, "g") for mol in ["H2", "CO", "CO2", "CH4", "C2H4", "O2"]]
    + [
        (mol, "l")
        for mol in [
            "H2O",
            "HCOOH",
            "CH3OH",
            "CH3COOH",
            "CH3CH2OH",
            "CH3CH2CH2OH",
            "CH3CHO",
            "C2H2O4",
        ]
    ]
    + [("e-", "")]
)

standard_redox = {
    "H": +1,
    "O": -2,
    "N": -3,
    "F": -1,
    "Cl": -1,
    "Br": -1,
    "I": -1,
    "S": -2,
    "Li": +1,
    "Na": +1,
    "K": +1,
    "Rb": +1,
    "Cs": +1,
    "e": -1,  # electrons come with a negative charge
    "-": +1,
    "+": -1  # If charges are counted this way, then
    # the sum of redox states in any compound is zero, which is convenient
}

dsH0 = {  # enthalpy of solvation / [kJ/mol], for T dependence of kH
    "ethanol": -19.5,  # temperarily taken value for CO2 below to check that the function runs properly! Couldn't find it for ethanol
    "CO2": -19.5,  # Carroll1991
}  # solvation enthalpy at 25C

kH0 = {  # Henry's Law constant of volatility in bar*l/mol
    "N2": 1660.87,
    "CO2": 29.87,
    "H2O": 0.0005791,
    "CH3CH2OH": 0.0047595,
    "Cl2": 10.411,
    "CO": 1078.33,
    "C2H4": 213.187,
    "C2H6": 52.553,
    "CH3CHO": 0.0714,
    "Ar": 728.80,
    "O2": 768.464,
    "CH4": 713.928,
    "He": 2726.81,
    "H2": 1289.037,
    "CH3OH": 0.00455,
    "CH3CH2CHO": 0.0769,
    "CH3CH2CH2OH": 0.00667,
    "HCOOH": 1.78e-4,
}  # all from Sander1999, I think. NIST WebBook quotes these values

aka = {  # dictionary of pseudonyms
    "C3H8O": "CH3CH2CH2OH",
    "CH3CH2CH2OH(l)": "CH3CH2CH2OH(aq)",
    "CH3CHO(l)": "CH3CHO(aq)",
    # ^ the standard states are liquid, but I have data for the aqueous
    "CO1": "CO",
}


def nu_to_rxn(nu, arrow="-->"):
    """
    --- arguments
    nu is a dictionairy containing the stoichiometric coefficients, such as
    {'CO2':-6, 'H2O':-6, 'C6H12O6':1, 'O2':6}
    --- return
    rxn is a string describing a reaction, such as
    '6 CO2 + 6 H2O -> C6H12O6 + 6 O2'
    """
    rxn = arrow
    for mol, n in nu.items():
        if type(n) is not int:
            n = np.round(n, decimals=2)
        if n > 0:
            if not rxn[-len(arrow) :] == arrow:
                rxn = rxn + " +"
            if n == 1:
                rxn = rxn + " " + mol
            else:
                rxn = rxn + " " + str(n) + " " + mol
        elif n < 0:
            if not rxn[0 : len(arrow)] == arrow:
                rxn = "+ " + rxn
            if n == -1:
                rxn = mol + " " + rxn
            else:
                rxn = str(-n) + " " + mol + " " + rxn
    return rxn


def rxn_to_nu(rxn, arrow=None):
    """
    --- arguments
    rxn is a string describing a reaction, such as
    '6 CO2 + 6 H2O -> C6H12O6 + 6 O2
    --- return
    nu is a dictionairy containing the stoichiometric coefficients, such as
    {'CO2':-6, 'H2O':-6, 'C6H12O6':1, 'O2':6}
    """

    nu = {}
    parts = rxn.split()
    if arrow is None:
        arrows = ["->", "-->"]
    elif type(arrow) is str:
        arrows = [arrow]
    lr = -1  # -1 for left of arrow, +1 for right of arrow
    n = 1  # nu[part]
    for part in parts:
        if part in arrows:
            lr = 1
            continue
        if part == "+":
            continue
        try:
            n = int(part)
            continue
        except ValueError:
            try:
                n = float(part)
                continue
            except ValueError:
                pass
        nu[part] = n * lr
        n = 1
    return nu


def get_formation_reaction(comp, out="nu", verbose=False):
    """
    Returns the formation reaction of comp, either as a string (out='string')
    or as a dictionairy of stoichiometric coefficents (out='nu')
    """
    nu = {comp: 1}
    elements = get_elements(comp)
    for atom, n in elements.items():
        ss = pure_states[atom]
        n_atom = get_elements(ss)[atom]
        if ss in nu:
            nu[ss] += -n / n_atom
        else:
            nu[ss] = -n / n_atom
    rxn = nu_to_rxn(nu)
    if verbose:
        print(rxn)
    if out == "string":
        return rxn
    return nu


def read_state(comp, verbose=False):
    match_state = re.search("\([a-z]+\)\Z", comp)
    if match_state is None:
        if verbose:
            print("can't read state for " + comp)
        return comp, None
    c = comp[: match_state.start()]
    s = match_state.group()[1:-1]
    return (c, s)


def get_cs(c, s):
    """
    cs is a compound with its state in parentheses. This function generates

    """
    if s is None or s == "":
        return c
    s = s.strip()
    c = c.strip()
    if not (s[0] == "(" and s[-1] == ")"):
        s = "(" + s + ")"
    cs = c + s
    return cs


def get_standard_state(comp, states={}, verbose=True, out="cs"):
    c, s = read_state(comp)
    if s is not None:
        if verbose:
            print(comp + " already has a state! Using state = " + s)
        ss = s
        comp = c
    else:
        try:
            ss = states[comp]  # compound with state
        except KeyError:
            if verbose:
                print(
                    "no state given for "
                    + comp
                    + ".\n"
                    + "Input it as, e.g. states={'"
                    + comp
                    + "':'aq'}"
                )
            if "+" in comp or "-" in comp:
                ss = "aq"
                if verbose:
                    print("I'll assume you meant " + ss)
            else:
                try:
                    ss = standard_states[comp]
                    if verbose:
                        print("Using its standard state, " + ss)
                except KeyError:
                    ss = None
    cs = get_cs(comp, ss)
    # print(cs)
    if out == "cs":
        return cs
    return comp, ss


def get_dfS(comp, T=None, verbose=True):
    """
    Get the change in entropy in the standard formation reaction of a compound
    """
    if T is not None:
        print("T-dependence of entropy not implemented. Using S0.")
    nu = get_formation_reaction(comp)
    dfS0 = 0
    for c, n in nu.items():
        #        print('c = ' + c + ', n = ' + str(n)) #for debugging
        try:
            dfS0 += n * S0[c]
        except KeyError:
            try:
                cs = get_standard_state(c, verbose=verbose)
                dfS0 += n * S0[cs]
            except KeyError:
                if verbose:
                    print(
                        "no standard entropy available for "
                        + c
                        + ". Couldn't get dfS for "
                        + comp
                        + ". Returning None."
                    )
                return None
    return dfS0


def get_dfH(comp, T=None, verbose=True):
    """
    Get the change in enthalpy in the standard formation reaction of a compound
    """
    if T is not None:
        print("T-dependence of enthalpy not implemented. Using H0.")
    try:
        dfH = dfH0[comp]
    except KeyError:
        try:
            cs = get_standard_state(comp, verbose=verbose)
            dfH = dfH0[cs]
        except KeyError:
            if verbose:
                print(
                    "no standard enthalpy available for " + comp + ". Returning None."
                )
            return None
    return dfH


def get_dfG(
    comp,
    T=None,
    dfG={},
    states={},
    T0=298.15,
    trycs=True,
    tryaka=True,
    verbose=True,
    vverbose=False,
):
    """
    Returns formation free energy of a compound in the specified state in kJ/mol

    e.g. get_dfG('H2O(l)') returns the standard free energy change of the
    reaction 'H2(l) + 1/2 O2(l) --> H2O(l)', which is dfG = -237.2 kJ/mol.

    This is a stubborn and robust function that really really wants to give
    you a free energy of formation for your compound.

    This compound first checks in the input dfG (enables it's automated use
    even when there's untabulated stuff in play.)
    Then it checks the tabulated value in the dictionairy dfG0 at the top
    of this module.
    Then it checks if it can calculate dfG from dfS and dfH.
    Then it checks if it can calcualte dfG from another state using e.g.,
    the henry's-law constant to calculate the free energy of acqueous species
    from the gas-phase free energy.
    Finally, it checks if the compound was input without its state, in which
    case it guesses a state and starts over.
    """
    dfG1 = dfG0.copy()
    dfG1.update(dfG)

    if T is not None:
        print("Temperature-dependent free energy not fully implemented.")
        try:
            dfH = get_dfH(comp, T=T)
            dfS = get_dfS(comp, T=T)
            dfG = dfH - T * dfS * 1e-3
            return dfG
        except KeyError:
            print("Something missing. returning None for dfG(" + comp + ").")
            return None

    if comp in dfG1:
        if verbose:
            print("found dfG(" + comp + ")")
        return dfG1[comp]

    if verbose:
        print("couldn't find dfG(" + comp + ")")

    # Try and get it from enthalpy and entropy:
    dfS = get_dfS(comp, verbose=vverbose)
    dfH = get_dfH(comp, verbose=vverbose)
    if dfS is not None and dfH is not None:
        dfGc = dfH - T0 * dfS * 1e-3
        print("returning dfH(" + comp + ") - T0 * dfS(" + comp + ")")
        return dfGc

    # Try and get it from other states of the same substance:
    c, s = read_state(comp)
    if s == "aq":
        cs_g = get_cs(c=c, s="g")
        dfGc_g = get_dfG(
            cs_g,
            T=T,
            dfG=dfG,
            states=states,
            T0=T0,
            trycs=trycs,
            tryaka=tryaka,
            verbose=verbose,
            vverbose=vverbose,
        )
        if dfGc_g is not None:  # then use the Henry's-law constant!
            kH = get_kH(c, tryG=False)  # tryG = False prevents the
            # infinite recursion that results  if it tries to get kH from dfG, that
            if kH is None:
                print("just using dfG(" + cs_g + ")")
                return dfGc_g
            dfGc = dfGc_g + R * T0 * 1e-3 * np.log(kH)
            if verbose:
                print("returning dfG(" + cs_g + ") + RTln(kH)")
            return dfGc

    print("couldn't get dfG for " + comp + ".")

    if tryaka:
        if comp in aka:
            print(comp + " is also known as " + aka[comp])
            return get_dfG(
                aka[comp],
                T=T,
                dfG=dfG,
                states=states,
                T0=T0,
                trycs=trycs,
                tryaka=False,
                verbose=verbose,
                vverbose=vverbose,
            )

    if trycs and s is None:  # likely, the input just forgot the state
        c, s = get_standard_state(comp, states=states, out="both", verbose=vverbose)
        cs = get_cs(c, s)
        return get_dfG(
            cs,
            T=T,
            dfG=dfG,
            states=states,
            T0=T0,
            trycs=False,
            tryaka=tryaka,
            verbose=verbose,
            vverbose=vverbose,
        )

    print("Returning None.")
    return None


def get_drG(nu, states={}, dfG={}, verbose=True):
    """
    returns standard free energy of formation in [kJ/mol] for a reaction given
    as a dictionary of stoichiometric coefficients or as a reaction string.
    """
    drG = 0
    if type(nu) is str:
        nu = rxn_to_nu(nu)
    for comp, n in nu.items():
        dfGc = get_dfG(comp, states=states, dfG=dfG, verbose=verbose)
        if dfGc is None:
            print(
                "no free energy of formation available for "
                + comp
                + ".\n"
                + "Input it as, e.g. dfG={'"
                + comp
                + "':0}.\n"
                + "for now, I'll assume it's zero"
            )
        else:
            drG += n * dfGc
    return drG


def p_vap(mol="H2O", T=298.15, unit="Pa"):
    """
    Returns the vapor pressure of a molecule at a given temperature, based on
    data in dfH0 and S0 dictionaries.
    """
    dH = (dfH0[mol + "(g)"] - dfH0[mol + "(l)"]) * 1e3
    dS = S0[mol + "(g)"] - S0[mol + "(l)"]

    if unit == "Pa":
        p0 = 1e5
    elif unit == "bar":
        p0 = 1
    elif unit == "mbar":
        p0 = 1000

    p = p0 * np.exp(-dH / (R * T) + dS / R)

    return p


def get_kH(
    comp=None, T=None, dfG={}, kH_0=None, dsH=None, T0=298.15, verbose=True, tryG=True
):
    """
    Returns the henry's-law vapor pressure of a molecule at a given temperature.
    Yet to be fully implemented.
    """
    if T is not None:
        if kH_0 is None:
            kH_0 = kH0[comp]
        if dsH is None:
            dsH = dsH0[comp]
        kH_T = kH_0 * np.exp(-dsH / R * (1 / T0 - 1 / T))
        return kH_T
    if comp in kH0:
        return kH0[comp]
    c, s = read_state(comp)
    if c in kH0:
        return kH0[c]
    if tryG:
        cs_aq = get_cs(c, "aq")
        dfG_aq = get_dfG(cs_aq, T=T, dfG=dfG, verbose=verbose)
        cs_g = get_cs(c, "g")
        dfG_g = get_dfG(cs_g, T=T, dfG=dfG, verbose=verbose)

        deG = dfG_g - dfG_aq  # free energy change of evaporation, i.e., aq --> g

        kH = np.log(deG * 1e3 / (R * T))
        print("got kH from dfG(" + cs_aq + ")")
        return kH
    print("couldn't get kH")


def get_oxidation_state(comp, atom="C", redox_states={}):
    """
    Gets the average oxidation state of a specified atom (default is carbon)
    in a chemical formula, with all the other elements having standard oxidation
    states or oxidation states specified here in redox_states.
    """
    redox = standard_redox.copy()
    redox.update(redox_states)
    if type(comp) is str:
        elements = get_elements(comp)
    else:
        elements = comp
    ro = 0  # the combined oxidation state of everything other than atom
    for element, n in elements.items():
        if element == atom:
            n_comp = n
            continue
        try:
            ro += n * redox[element]
        except KeyError:
            print(
                "I don't know the oxidation state of "
                + element
                + "\n"
                + "Input it in as redox_states = {atom: oxidation_state}.\n"
                + "For now, I'll assume it's zero."
            )

    return -ro / n_comp


def get_rxn_EC(
    product="CH4",
    reactant="CO2",
    atom="C",
    redox_states={},
    extra=["H+", "H2O"],
    out="string",
):
    """
    Generates a balanced electrochemical reaction
    inputs
        - product : formula for product compound
        - reactant : formula for reactant compound
        - atom : we assume only this element changes oxidation state
        - redox_states : dict for which each key is a spectator atom who's
                        redox state is assumed to be the corresponding value.
                        Default is at the top of this script
        - extra : will be the compounds used to balance the spectator atoms.
                    NOT YET IMPLEMENTED. H2O and H+ are always used now.
        - out : specifies output. default is reaction string.
                out = 'nu' for stoichiometric coefficients of reaction
    """
    redox = standard_redox.copy()
    redox.update(redox_states)

    elements_r = get_elements(reactant)
    natoms_r = elements_r[atom]
    ro_r = get_oxidation_state(elements_r, atom=atom, redox_states=redox_states)
    elements_p = get_elements(product)
    natoms_p = elements_p[atom]
    ro_p = get_oxidation_state(elements_p, atom=atom, redox_states=redox_states)

    nu = {}  # balanced reaction will go here
    if product == reactant:
        print("product equals reactant. No reaction.")
        return nu

    n_atoms = natoms_r * natoms_p / gcd(natoms_r, natoms_p)
    # total number of specified atom in balanced reaction

    nu[reactant] = int(-n_atoms / natoms_r)
    nu[product] = int(n_atoms / natoms_p)

    nu["e-"] = n_atoms * (ro_p - ro_r)

    # this can be done in a cool general way, but, assuming I've just got H, O, and atom:
    if "O" in get_elements(nu):
        if "H2O" in nu:
            nu["H2O"] += -get_elements(nu)["O"]
        else:
            nu["H2O"] = -get_elements(nu)["O"]
    if "H" in get_elements(nu):
        if "H+" in nu:
            nu["H+"] += -get_elements(nu)["H"]
        else:
            nu["H+"] = -get_elements(nu)["H"]

    if not set(get_elements(nu).keys()).issubset({atom, "O", "H", "+", "-", "e"}):
        print(
            "function Chem.get_rxn_EC: WARNING!\n"
            + "Some of the atoms in your reaction aren't implemented.\n"
            + "The reaction won't be balanced."
        )

    rxn = nu_to_rxn(nu)

    print(rxn)
    if out == "string":
        return rxn
    return nu


def get_rxn_cell(rxn_an, rxn_cat, out="string"):
    if type(rxn_an) is str:
        rxn_an = rxn_to_nu(rxn_an)
    if type(rxn_cat) is str:
        rxn_cat = rxn_to_nu(rxn_cat)
    n_an = int(rxn_an["e-"])  # positive for forward reaction
    n_cat = int(rxn_cat["e-"])  # negative for forward reaction
    n = (
        -n_an * n_cat / gcd(n_an, n_cat)
    )  # positive if both are forward or both backwards
    nu = {}
    for comp, i in rxn_an.items():
        if comp in nu:
            nu[comp] += n / n_an * i  # positive if rxn_cat is forwards
        else:
            nu[comp] = n / n_an * i
    for comp, i in rxn_cat.items():
        if comp in nu:
            nu[comp] += -n / n_cat * i  # positive if rxn_an is forwards
        else:
            nu[comp] = -n / n_cat * i
    rxn = nu_to_rxn(nu)
    print(rxn)
    if out == "string":
        return rxn
    return nu


def get_rxn_c(
    reactant="CH4",
    product="CO2",
    atom="C",
    redox_states={},
    extra=["H+", "H2O"],
    out="string",
):
    """
    Returns the combustion reaction for a product.
    """

    nu_cat = get_rxn_EC(product="H2O", reactant="O2", atom="O", out="nu")
    nu_an = get_rxn_EC(
        product=product,
        reactant=reactant,
        redox_states=redox_states,
        extra=extra,
        out="nu",
    )
    nu_c = get_rxn_cell(nu_an, nu_cat, out="nu")
    if out == "string":
        return nu_to_rxn(nu_c)
    else:
        return nu_c


def get_standard_potential(nu, pH=0, T=298.15, states={}, dfG={}, verbose=True):
    """
    returns equilibrium potential in [V] on SHE scale for a reaction given
    as a dictionary of stoichiometric coefficients or as a rxn string.
    If nu['H+'] == nu['e-'], then leave pH=0 to get the pH-independent
    standard potential on the RHE scale
    """
    if type(nu) is str:
        nu = rxn_to_nu(nu)
    drG = get_drG(nu, states=states, dfG=dfG, verbose=verbose)
    n_el = nu["e-"]
    E0 = drG * 1e3 / (n_el * Far)
    if "H+" in nu.keys():
        n_H = nu["H+"]
        E0 = E0 - n_H / n_el * np.log(10) * R * T / Far * pH
    return E0


def get_dcG(
    reactant="CH4",
    product="CO2",
    atom="C",
    redox_states={},
    extra=["H+", "H2O"],
    states={},
    dfG={},
    verbose=True,
):
    nu = get_rxn_c(
        product=product,
        reactant=reactant,
        redox_states=redox_states,
        extra=extra,
        out="nu",
    )
    dcG = get_drG(nu, states=states, dfG=dfG, verbose=verbose)
    if not nu[reactant] == -1:
        dcG = -dcG / nu[reactant]
    return dcG
