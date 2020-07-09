# %%


from matplotlib import pyplot as plt
from EC_MS import Molecule

from EC_MS.utils.extraction_class import Extraction

# %%
Molecule("ethanol").plot_spectrum()

# %%

plt.show()

print(Extraction.__str__)
