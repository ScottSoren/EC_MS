# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 19:08:46 2020

@author: scott
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from EC_MS import Spectrum

plt.close("all")

"""
This script demonstrates the Spectrum class. All its methods are quite simple,
except for the initialization method whcih parses the file. This script
shows you what the methods do. There are usually more than one way of doing 
things, so I sometimes show both, and you can choose what's most intuitive for you.
"""

folder = os.path.abspath(".")  # for this demo, I just copied the data files into here.
file = folder + os.sep + "1_SEM_5E-6.dat"

# the data is naturally imported as a Spectrum object.
spectrum = Spectrum(file, data_type="PVMS")
# ^ it complains because it can't find the timestamp in the file name. You can just ignore that.
print("\n" + "-" * 20 + "\n")
# the spectrum knows its file and folder.
print("spectrum.file = " + spectrum.file + "\nspectrum.folder = " + spectrum.folder)

ax = spectrum.plot(color="k")  # kwargs like 'color' are fed to plt.plot()
ax.set_yscale("log")

# Manually get max and integral integrate:
x, y = spectrum.get_signal(
    Mspan=[27.5, 28.5]
)  # spectrum.get_signal selects the part of the spectrum in Mspan
int_M28_naive = np.trapz(y, x)
print(
    "The integrated intensity of the primary N2 peak is " + str(int_M28_naive) + " C."
)
print("... maximum = " + str(max(y)))
# but this is not completely right as it ignores the background, as seen here:
ax.fill_between(x, y, np.zeros(x.shape), color="r")

# the Spectrum.get_integral() function does so automatically. It plots if given an axis
int_M28 = spectrum.get_integral(Mspan=[27.5, 28.5], M_bg=[70, 80], ax=ax, color="b")
# ^ you can see from the blue fill where it set the background.
print(
    "The background-corrected integrated intensity of the primary N2 peak is "
    + str(int_M28)
    + " C."
)


# we can also do a global background correction:
spectrum.set_background(M_bg=[70, 80])
# now the max is a tiny bit lower. We can also automatically get the max, like this:
y_max = spectrum.get_max(Mspan=[27.5, 28.5])
print("... background-corrected maximum = " + str(y_max) + " A.")

# we can plot it as before, or manually:
x = spectrum.x
y = spectrum.y
fig, ax = plt.subplots()
ax.plot(x, y, "k")
ax.set_yscale("log")

# ah, but it looks like shit because we're plotting it on a log scale. Luckily, we can undo the background subtraction!
spectrum.reset()  # undoes background subtraction
spectrum.plot(ax=ax, Mspan=[10, 80], color="b", linestyle="--")

plt.show()
# that concludes this little demo.
