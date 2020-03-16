# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:34:09 2020

@author: scott
"""

import numpy as np
from matplotlib import pyplot as plt
from EC_MS import Chip, Chem

chip = Chip()

print("design length = " + str(chip.l_cap * 1e3) + " mm")


T0 = 273.15
T_vec = np.linspace(0, 100, 101) + T0


fig, ax = plt.subplots()
ax.set_xlabel("Temperature / [deg C]")
ax.set_ylabel("capilarry flux / [nmol/s]")

p_colors = [(1e5, "k"), (2e5, "b"), (2.5e5, "g"), (3e5, "r")]


for p, color in p_colors:

    n_dot_vec = np.array([])
    for T in T_vec:
        # chip.T = T
        # chip.p = p
        n_dot = chip.capillary_flow("air", p=p, T=T) / Chem.NA
        n_dot_vec = np.append(n_dot_vec, n_dot)

    ax.plot(T_vec - T0, n_dot_vec * 1e9, color=color, label=str(p))

ax.legend()
fig.savefig("air.png")
