# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:50:10 2016
Most recently edited: 16J27

@author: scott

This module will include functions for modelling mass transport
and fitting time response in the EC-MS setup.

See Scott's MSc thesis, chapter 2 and section 3.3, for discussion inc. prior
implementation in Matlab


"""


from __future__ import division, print_function

import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint

from . import Chem
from .Molecules import Molecule
from .Plotting import plot_operation


def fit_exponential(t, y):
    """
    A previous attempt at this had used scipy.optimize.minimize.
    """
    t = t - t[0]  # zero time axis
    tau_i = t[-1] / 10  # guess at time constant
    # tau_i = t[-1]      #often can't solve with this guess. A smaller tau helps.
    y0_i = y[-1]  # guess at approach value
    y1_i = y[0]  # guess at true initial value
    pars_i = [tau_i, y0_i, y1_i]

    def exp_fun(x, tau, y0, y1):
        z = y0 + (y1 - y0) * np.exp(-x / tau)
        #        print([tau,y0,y1]) #for diagnosing curve_fit problems
        return z

    pars, pcov = curve_fit(exp_fun, t, y, p0=pars_i)
    #    pars = [tau, y0, y1]

    return pars


def fit_step(t, y, tpulse=0, step="fall", ax=None, spec="r--", label=None, verbose=1):
    """
    use step='rise' to fit onset and step='fall' to fit tail
    assumes that data starts from the start of the pulse
    """
    if verbose:
        print("\n\nfunction 'fit_step' at your service!\n")

    # zero time axis
    t0 = t[0]
    t = t - t0
    print("t0 = " + str(t0))
    if type(t) is list:
        t = np.array(t)  # 17B02
    if step == "fall":
        I_tail = np.array([I for (I, t_I) in enumerate(t) if tpulse < t_I])
        #       print(I_tail)
        t_fit = t[I_tail] - tpulse
    elif step == "rise":
        if tpulse == 0:
            tpulse = t[-1]
        I_tail = np.array([I for (I, t_I) in enumerate(t) if t_I < tpulse])
        t_fit = t[I_tail]
    else:
        print("use step='rise' to fit onset and step='fall' to fit tail")

    pars = fit_exponential(t_fit, y[I_tail])

    if ax:
        tau = pars[0]
        y0 = pars[1]
        y1 = pars[2]
        y_fit = y0 + (y1 - y0) * np.exp(-t_fit / tau)
        t_fit = t_fit + t0  # put time axis back
        if step == "fall":
            t_fit = t_fit + tpulse
        ax.plot(t_fit, y_fit, spec, label=label)
        if label:
            if label == "tau":
                label = "tau = {0:5.2f} s".format(tau)
            I_text = int(len(t_fit) / 2)
            t_text = t_fit[I_text]
            y_text = y_fit[I_text]
            ax.text(t_text, y_text, label, color=spec[0])

    if verbose:
        print("tau = " + str(pars[0]) + " s")
        print("\nfunction 'fit_step' finished!\n\n")
    return pars


def stagnant_diffusion_ode(C, T, pars):  # Note that C comes before T here!
    """

    Scott's master p. 52 and appendix C.  Z is the new X.

    returns rate of change dC/dT of concentration profile for
    non-dimensionalized stagnant sniffer diffusion problem.
    C = C(X) where X goes from 0 (membrane) to 1 (electrode)
    T is time non-dimensionalized on the diffusion scale t0 = L²/D
    pars:
    [0] alpha = h*L/D is the system parameter.
    [1] J_fun returns the flux from the electrode as a unction of T. The flux
    scale J0 is used to define the concentration scale, C0 = J0*L/D

    #modified 17B02 to enable carrier gas introduction of element using Cg

    # pars as a list rather than a dict is slightly faster (see fig06.out),
    # but I like the dict so that I can remember what's what.
    """
    dZ = pars["dZ"]  # assigning this here replaces two lookups with one.

    C_N = C[-1] + pars["J_fun"](T) * dZ  # boundary condition dC/dZ = J(T) at electrode
    C_ = (
        C[0] - pars["alpha"] * (C[0] - pars["Cg"]) * dZ
    )  # boundary condition dC/dZ = -alpha*(C-Cg) at membrane
    C_up = np.append(C[1:], C_N)
    C_down = np.append(C_, C[:-1])
    d2CdZ2 = (C_up - 2 * C + C_down) * pars["1/dZ**2"]  # second derivative of C wrt Z
    dCdT = d2CdZ2  # Fick's second law

    return dCdT


def solve_stagnant(
    alpha=1, J_fun=None, Cg=1, Tspan=[0, 10], startstate="zero", N=30, flux=1, verbose=1
):
    """solves the stagnant sniffer partial differential equations.
    pars[0][0] is alpha = h*L/D is the system parameter.
    pars[0][1] is J_fun. Returns the flux from the electrode as a unction of T.
    Tspan is [Tstart, Tfinish] on the diffusion timescale t0 = L²/D
    C0 = C0(X) is start concentration profile. If size(C0) = 1, assumes a
    uniform concentration.
    N is descretization (read from C0)
    flux = 0 to return entire concentration profile (on C0 = J0*L/D scale)
    flux = 1 to return flux through membrane (on J0 scale)
    """
    if verbose:
        print("\n\nfunction 'solve_stagnant' at your service!\n")

    if startstate == "zero":
        C0 = np.zeros([N])
    elif startstate == "steady":
        C0 = 1 / alpha + np.linspace(0, N) / (N + 1)  # Scott's MSc thesis, p. 53
    elif startstate == "saturated":
        C0 = np.ones([N]) * Cg
    elif np.size(startstate) == 1:
        C0 = np.ones([N]) * startstate
    else:
        C0 = startstate
        N = np.size()

    if np.size(Tspan) == 2:
        Tspan = np.linspace(Tspan[0], Tspan[1], 100)

    dZ = 1 / N  # N+1? N-1? I can never remember what's most 'correct'

    pars = {"alpha": alpha, "J_fun": J_fun, "Cg": Cg, "dZ": dZ, "1/dZ**2": 1 / dZ ** 2}

    CC = odeint(stagnant_diffusion_ode, C0, Tspan, args=(pars,))
    # 16J18_02h10: this crashes the kernel, I don't know why... 18h56 found it! c before t!

    J = (
        CC[:, 1] - CC[:, 0]
    ) * N  # (positive) J = dC/dX with dC = C0 - C_ and dZ = 1 / N
    # J is a function of T

    if verbose:
        print("solution shape: " + str(np.shape(CC)))
        print("\nfunction 'solve_stagnant' finished!\n\n")
    if flux:
        return Tspan, J
    else:
        return Tspan, CC


def stagnant_pulse(*args, **kwargs):
    print(
        "\n\n'stagnant_pulse' has been renamed 'stagnant_operator'. Remember that next time!"
    )
    return stagnant_operator(*args, **kwargs)


def stagnant_operator(
    tj=None,
    tpulse=10,
    tspan=None,
    j_el=-1,
    L=100e-6,
    A=0.196e-4,
    q0=1.5e15 / Chem.NA,
    p_m=1e5,
    mol="H2",
    p_gas=0,
    normalize=False,
    D=None,
    kH=None,
    n_el=None,
    Temp=None,
    unit="pmol/s",
    flux_direction="out",
    verbose=True,
    ax=None,
    plot_type=None,
    startstate="zero",
    N=30,
    colormap="plasma",
    aspect="auto",
):
    """
    Models a pulse of current towards a specified product in our EC-MS setup.
    Theory in chapter 2 of Scott's masters thesis.
    all arguments are in pure SI units. The electrode output can either be given
    as a steady-state square pulse of electrical current (tpulse, j_el, n_el),
    or as a measured current (tj[1]) as a function of time (tj[0])

    #tj[1] should have units A/m^2. 1 mA/cm^2 is 10 A/m^2

    #17B02: p_gas is the partial pressure of the analyte in the carrier gas.
    # this enables, e.g., CO depletion modelling.
    """
    if verbose:
        print("\n\nfunction 'stagnant_operator' at your service!\n")
    if type(mol) is str:
        mol = Molecule(mol)
    if Temp is not None:
        mol.set_temperature(Temp)
    else:
        Temp = 298.15  # standard temperature in K
    if D is None:
        D = mol.D
    if kH is None:
        kH = mol.kH
    if n_el is None and not normalize:
        n_el = mol.n_el
    if tspan is None:
        if tj is None:
            tspan = [-0.1 * tpulse, 1.2 * tpulse]
        else:
            tspan = [tj[0][0], tj[0][-1]]

    h = kH * Chem.R * Temp * q0 / (p_m * A)  # mass transfer coefficeint

    alpha = L * h / D  # system parameter

    # non-dimensional scales:
    t0 = L ** 2 / D
    if tj is None:
        if normalize:
            j0 = 1
        else:
            j0 = j_el / (n_el * Chem.Far)
    else:
        t = tj[0]
        if normalize:
            j0 = 1
            j = tj[1] / np.max(np.abs(tj[1]))
        else:
            j = tj[1] / (n_el * Chem.Far)  # A/m^2 --> mol/(m^2*s)
            j0 = max(np.abs(j))
    c0 = j0 * L / D
    tau = L ** 2 / (2 * D) + L / h
    # from the approximate analytical solution, Scott's thesis appendix D

    Tpulse = tpulse / t0
    Tspan = (
        np.linspace(tspan[0], tspan[1], 1000) / t0
    )  # why do I give so many time points?

    if tj is None:

        def J_fun(T):
            if T < 0:
                return 0
            if T < Tpulse:
                return 1
            return 0

    else:
        T_in = t / t0
        J_in = j / max(np.abs(j))
        # print('max(J_in) = ' + str(max(J_in)))
        def J_fun(T):
            if T < T_in[0]:  # assume no current outside of the input tj data
                return 0
            if T < T_in[-1]:
                return np.interp(T, T_in, J_in)
            return 0

    c_gas = p_gas / (Chem.R * Temp)
    cg = c_gas / kH  # concentration analyte in equilibrium with carrier gas, 17B02
    Cg = (
        cg / c0
    )  # non-dimensionalized concentration analyte at equilibrium with carrier gas, 17B02

    # pars = ([alpha, J_fun, Cg],) #odeint needs the tuple. 17A12: Why ?!
    [T, CC] = solve_stagnant(
        alpha=alpha,
        J_fun=J_fun,
        Cg=Cg,
        Tspan=Tspan,
        startstate=startstate,
        flux=False,
        N=N,
    )
    cc = CC * c0
    t = T * t0
    j = h * (cc[:, 0] - cg)  # mass transport at the membrane
    # j1 = D * (cc[:,1] - cc[:,0])
    # fick's first law at the membrane gives the same j :)
    if verbose:
        print(
            "q0 = "
            + str(q0)
            + " mol/s, h = "
            + str(h)
            + " m/s, alpha = "
            + str(alpha)
            + ", j0 = "
            + str(j0)
            + " mol/(m^2*s), max(j)/j0 = "
            + str(max(j) / j0)
            + ",  t0 = "
            + str(t0)
            + " s, c0 = "
            + str(c0)
            + " mol/m^3"
            + ", tau (analytical) = "
            + str(tau)
            + " s"
            + ", cg = "
            + str(cg)
            + " mM"
        )

    # get ready to plot:
    N = np.shape(cc)[1]
    z = np.arange(N) / (N - 1) * L
    # this will only be used for heatmap, so it's okay if dx isn't quite right.
    if "cm^2" not in unit:
        j = j * A
    if unit[0] == "u":
        j = j * 1e6
    elif unit[0] == "n":
        j = j * 1e9
    elif unit[0] == "p":
        j = j * 1e12
    if flux_direction == "in":
        j = -j

    if normalize:
        s_int = np.trapz(j, t)
        if verbose:
            print("normalizing from area = " + str(s_int))
        j = j / s_int

        # plotting was moved on 17G30 some legacy code here:
    if plot_type is not None and ax is not None:
        print("We recommend you plot seperately, using the function 'plot_operation'.")
        axes = plot_operation(
            cc=cc,
            t=t,
            z=z,
            j=j,
            ax=ax,
            plot_type=plot_type,
            colormap=colormap,
            aspect=aspect,
            verbose=verbose,
        )
        if verbose:
            print("\nfunction 'stagnant_operator' finished!\n\n")
        return t, j, axes

    results = {"t": t, "z": z, "j": j, "cc": cc, "dimensions": "tz"}

    return results


def flow_diffusion_ode(C, X, pars):
    """
    Scott's master, p. 60. X is the new Y and Z is the new X.
    """
    C_N = C[-1]
    C_ = C[0] - pars["alpha"] * (C[0] - pars["Cg"]) * pars["dZ"]
    C_up = np.append(C[1:], C_N)
    C_down = np.append(C_, C[:-1])
    d2CdZ2 = (C_up - 2 * C + C_down) * pars["1/dZ**2"]
    # I assume multiplication to be faster than division
    dCdX = d2CdZ2 * pars["1/beta"]
    return dCdX


def solve_flow(
    alpha=1, beta=1, Cg=0, N=30, flux=False, verbose=True, Xspan=[0, 1], C0="uniform"
):
    """
    This solves the flow ODE and returns either flux through membrane as a
    function of position (Xspan and J),
    or the concentration profile (Xspan and CC)

    It assumes steady state. I think I can use this and then convolute if I
    need time dependence.
    """

    if verbose:
        print("\n\nfunction 'solve_flow' at your service!\n")

    if C0 == "uniform":
        C0 = np.array([1] * N)
        # nothing else really makes sense, since c0 defines a scale.

    if np.size(Xspan) == 2:
        Xspan = np.linspace(Xspan[0], Xspan[1], 100)

    dZ = 1 / N  # N+1? N-1? I can never remember what's most 'correct'

    pars = {
        "alpha": alpha,
        "1/beta": 1 / beta,
        "Cg": Cg,
        "dZ": dZ,
        "1/dZ**2": 1 / dZ ** 2,
    }

    CC = odeint(flow_diffusion_ode, C0, Xspan, args=(pars,))
    # 16J18_02h10: this crashes the kernel, I don't know why... 18h56 found it! c before t!

    J = (
        CC[:, 1] - CC[:, 0]
    ) * N  # (positive) J = dC/dZ with dC = C0 - C_ and dX = 1 / N
    # J is a function of X.

    if verbose:
        print("solution shape: " + str(np.shape(CC)))
        print("\nfunction 'solve_flow' finished!\n\n")
    if flux:
        return Xspan, J
    else:
        return Xspan, CC
    pass


def flow_operator(
    mode="steady",  # in steady mode it's not really an operator.
    system="chip",  #
    A_el=0.196e-4,
    A=0.196e-4,
    q0=1.5e15 / Chem.NA,
    Temp=None,  # universal pars
    L=100e-6,
    w=5e-3,
    w2=5e-3,
    F=1e-9,  # geometry pars   #w and w2 changed from 0.5e-3 to 5e-3 on 17K28.
    c0=None,
    j0=None,
    j_el=-1,  # inlet flow pars
    p_m=1e5,  # chip pars
    phi=0.5,
    dp=20e-9,
    Lp=100e-6,  # DEMS pars
    mol="H2",
    D=None,
    kH=None,
    n_el=None,
    M=None,  # mol pars
    p_gas=0,
    normalize=False,
    unit="pmol/s",
    flux_direction="out",
    N=100,  # solver pars
    verbose=True,
):
    """
    Follows the recipe in Scott's MSc, page 61, for calculating collection
    efficiency in a flow system by solving a differential equation. This
    can be used to compare different types of EC-MS
    """
    if verbose:
        print("\n\nfunction 'flow_operator' at your service!\n")
    if type(mol) is str:
        mol = Molecule(mol)
    if Temp is not None:
        mol.set_temperature(Temp)
    else:
        Temp = 298.15  # standard temperature in K
    if D is None:
        D = mol.D  # diffusion constant in electrolyte / [m^2/s]
    if kH is None:
        kH = mol.kH  # dimensionless henry's law constant
    if M is None:
        M = Chem.Mass(mol.name) * 1e-3  # molar mass / [kg/mol]
    # print(c0)
    if n_el is None and c0 is None and not normalize:
        n_el = mol.n_el
    if c0 is None:
        if j0 is None:
            j0 = j_el * A_el / (n_el * Chem.Far)
        c0 = j0 / F
        # concentration is production rate over flow rate: [mol/s] / [m^3/s)] = [mol/m^3] )

    if system == "chip":
        h = kH * Chem.R * Temp * q0 / (p_m * A)  # Scott's MSc, page 50
    elif system == "DEMS":
        h = (
            kH * phi * dp / (3 * Lp) * np.sqrt(8 / np.pi * Chem.R * Temp / M)
        )  # Scott's MSc, page 49

    v0 = F / (L * w2)
    alpha = h * L / D
    beta = v0 * L ** 2 / (D * w)  # There is a mistake in Scott's MSc page60!
    # There I got the non-dimensionalization wrong, and wrote beta = v0*w**2/(D*L)
    # in fact, beta = v0*L**2/(D*w)

    Xspan = np.linspace(0, 1, 1000)
    X, CC = solve_flow(alpha=alpha, beta=beta, Xspan=Xspan, N=N, verbose=verbose)

    x = X * w
    cc = CC * c0
    j = cc[:, 0] * h  # mol/m^3 * m/s = mol/(m^2*s)
    Z = np.linspace(0, 1, N)
    z = Z * L

    eta_m = 1 - np.trapz(CC[-1, :], Z)  # Scott's MSc, page 61
    eta_m_check = (
        w2 * np.trapz(j, x) / (c0 * F)
    )  # m*mol/(m^2*s)*m / ((mol/m^3)*m^3/s) = 1

    qm = c0 * F * eta_m

    if verbose:
        print("portion not escaped = " + str(eta_m))
        print("portion collected = " + str(eta_m_check) + "\n\n")
    if system == "chip":
        eta_v = 1
    elif system == "DEMS":
        p_w = Chem.p_vap(mol="H2O", T=Temp, unit="Pa")
        M_w = Chem.Mass("H2O") * 1e-3
        j_w = (
            A
            * p_w
            / (Chem.R * Temp)
            * phi
            * dp
            / (3 * Lp)
            * np.sqrt(8 / np.pi * Chem.R * Temp / M_w)
        )
        eta_v = q0 / (qm + j_w)  # fixed 17H10

    eta = eta_m * eta_v

    if verbose:
        print(
            "q0 = "
            + str(q0)
            + " mol/s, h = "
            + str(h)
            + " m/s, alpha = "
            + str(alpha)
            + ", j0 = "
            + str(j0)
            + " mol/s, max(c)/c0 = "
            + str(np.max(np.max(cc)) / c0)
            + ",  kH = "
            + str(kH)
            + ", eta = "
            + str(eta)
            + ", mol = "
            + str(mol.name)
            + ", system = "
            + str(system)
            + ""
            + ", beta = "
            + str(beta)
            + ""
            + ", v0 = "
            + str(v0)
        )

    if verbose:
        print("\nfunction 'flow_operator' at finished!\n\n")

    results = {
        "x": x,
        "z": z,
        "j": j,
        "cc": cc,
        "eta_m": eta_m,
        "eta_v": eta_v,
        "eta": eta,
        "dimensions": "xz",
    }
    return results


def delta_response(
    L=100e-6,
    q0=1e15 / Chem.NA,
    mol="H2",
    D=None,
    kH=None,
    n_el=None,
    A=0.196e-4,
    p_m=1e5,
    Temp=298.15,
    verbose=True,
    tspan="auto",
    N_t=1000,
):
    """
    Returns the output when a stagnant_operator operates on a delta function.
    There's probably a much smarter way to do it, but for now I'll just do
    triangle pulse of width tau/250
    """
    if D is None or kH is None:
        if type(mol) is str:
            mol = Molecule(mol)
        if D is None:
            D = mol.D
        if kH is None:
            kH = mol.kH
    if verbose:
        print("calculating a delta function response.")
    h = kH * Chem.R * Temp * q0 / (p_m * A)  # mass transfer coefficeint
    tau = L / h + L ** 2 / (2 * D)
    if tspan == "auto":
        tspan = [0, 4 * tau]

    t = np.linspace(tspan[0], tspan[1], N_t)
    j = np.append(np.array([1]), np.zeros(N_t - 1))
    tj = [t, j]
    print(type(tj))
    return stagnant_pulse(
        tj=tj,
        normalize=True,
        tspan=tspan,
        L=L,
        A=A,
        q0=q0,
        p_m=p_m,
        D=D,
        kH=kH,
        n_el=n_el,
        Temp=Temp,
        verbose=True,
        plot_type=None,
    )


if __name__ == "__main__":

    pass

    ##
