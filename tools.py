
from src.tools.xls2npz import *
from src.ct.def_ct_tools import *
import matplotlib.pyplot as plt
from src.tools.def_DNS import *
from src.core.def_GPSA import *
from src.ck.def_cheminp import *


fld = '/Users/xianggao/GPS data/prj_oxy'
detailed_folder = os.path.join(fld, 'detailed','mech')
sk_folder = os.path.join(fld, 'ske19')
species_kept = 'CH3, CH4, CH3O, CH2O, HCO , CO, CO2, C2H6, C2H5, C2H4, C2H3, H, O, O2, H2O2, OH, H2, H2O, HO2'.replace(' ','').split(',')
skeletal(detailed_folder, sk_folder, species_kept)




