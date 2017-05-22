'''
Most recently edited 16I23. Imports the EC_MS package

For a description of this package, see the included NOTES__EC_MS.text


@author: Scott
'''


from EC_MS.Data_Importing import import_data, numerize
from EC_MS.Combining import *
from EC_MS.Plotting import *
from EC_MS.EC import *
from EC_MS.Integrate_Signals import integrate_pulses, step_vs_potential, integrate_transient
import EC_MS.Isotope_Ratio
import EC_MS.ML_strip
import EC_MS.Chem as Chem

