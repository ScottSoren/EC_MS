'''
Most recently edited 16I23. Imports the EC_MS package

For a description of this package, see the included NOTES__EC_MS.text


@author: Scott
'''

__version__ = '0.5.0.dev0'
__title__ = 'EC_MS'
__description__ = 'Analysis tools for electrochemistry and mass spectrometry and a lot in between'
__url__ = "https://github.com/ScottSoren/EC_MS"
__author__ = "Soren B. Scott"
__email__ = "sbs@spectroinlets.com"

__license = "MIT"

print('\n' + '-'*10 + '  Importing: EC_MS ' + __version__ + ' ' + '-'*10)
print('from ' + __file__ + '\n\n')


from EC_MS.Data_Importing import *
from EC_MS.Combining import *
from EC_MS.Plotting import *
from EC_MS.EC import *
from EC_MS.Quantification import *
from EC_MS.Molecules import *
from EC_MS.Chips import *
from EC_MS.Object_Files import *
from EC_MS.Integrate_Signals import *
from EC_MS.Time_Response import *
from EC_MS.Calibration import *
from EC_MS.Electrolytes import *
from EC_MS.Datapoints import *
from EC_MS.Potentiostat import *
from EC_MS.dataset import *
from EC_MS.PVMassSpec import *
from EC_MS.patches import *
import EC_MS.Chem as Chem
