# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 16:03:33 2020

@author: scott
"""


import numpy as np
from scipy.optimize import minimize
from scipy.integrate import odeint
from scipy.interpolate import interp1d


def carbonic_ode(S, t, pars):
    """
    Equations from Scott's Thesis page 53
    """
    k = pars[0]  # rate constant / s^-1
    alpha = pars[1]  # H2(16)O / H2(18)O ratio

    S44, S46, S48 = S[0], S[1], S[2]

    dS44 = k * (-2 / 3 * (1 - alpha) * S44 + 1 / 3 * alpha * S46)
    dS46 = k * (2 / 3 * (1 - alpha) * S44 - 1 / 3 * S46 + 2 / 3 * alpha * S48)
    dS48 = k * (1 / 3 * (1 - alpha) * S46 - 2 / 3 * alpha * S48)

    return [dS44, dS46, dS48]


def solve_carbonic_burst(k=0.026, alpha=0.27, tspan=[0, 60]):
    """
    Returns the partial concentrations at M44, M46, and M48 following an
    initial burst of CO(16) oxidation given:
        g = the H2(18)O/H2(16)O ratio
        k = the rate constant for H2O + CO2 --> H2CO3 in s^-1
    """
    print("k = " + str(k))
    print("alhpa = " + str(alpha))
    S0 = np.array([alpha, 1 - alpha, 0])
    pars = [k, alpha]
    if len(tspan) == 2:
        tspan = np.linspace(tspan[0], tspan[-1], 200)
    SS = odeint(carbonic_ode, S0, tspan, args=(pars,))
    return SS
