'''
@scott for EC_MS
created 16K02
based on theory in section 3.1 of Scott's MSc thesis
'''

import numpy as np

if os.path.split(os.getcwd())[1] == 'EC_MS':   
                       #then we're running from inside the package
    import Chem
    from Molecules import Molecule
else:
    from . import Chem
    from .Molecules import Molecule
#where should equilibrium constants be stored? In electrolyte objects?
#yeah... 

K_w = 1.0e-14

class Electrolyte(name):
    def __init__(self, name):
        self.name = name
        if name = 'phosphate':
            self.npK_a = np.array([2.15, 7.20, 12.32])
            self.z0 = 0
        elif name = 'carbonate':
            self.npK_a = np.array([3.6, 10.32])
            self.nK_eq = 1.7e-3   
            self.nK_h = 29.6  #bar/M
            self.z0 = 0
        self.K_a = np.power(10, -pKa)


def get_pH_sF(s, F, electrolyte='phosphate')
    '''
    here, s is solution strength (cation concentration)
    and F is the concentration of (anionic) species in 
    acid-base (ie, 'Fast') equilibruim
    '''
    ele = Electrolyte(electrolyte)
    K_a = ele.K_a

    def s_of_pH(pH):
        x = np.power(10, -pH)
        y = K_w/x
        n_base = range(len(K_a))
        q = sum(np.product(K_a[:i]) * np.power(x, i)  for i in n_base)
        anions = -sum((z0 - i)*np.product(K_a[:i] * np.power(x,i)   for i in range(n_base))/q
        s = y + anions - x
        return s

    x = np






