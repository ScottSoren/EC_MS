
'''
Imports Chem package, from within EC_MS package. Want to be able to use
Chem package both in scripts in EC_MS and in scripts that import EC_MS
there must be a better way. Need to ask!
'''

import os 

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the EC_MS package
    from Chem.MolarMasses import Mass
    from Chem.PhysCon import *
else:                           #then we use relative import
    from EC_MS.Chem.MolarMasses import Mass
    from EC_MS.Chem.PhysCon import *

