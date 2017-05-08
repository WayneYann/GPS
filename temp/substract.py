import os

from src.tools.xls2npz import *
from src.ct.def_ct_tools import *
import matplotlib.pyplot as plt
from src.tools.def_DNS import *
from src.core.def_GPSA import *
from src.post.def_painter import plot_GPedge, plot_GPedge_mf
from src.post.def_plt_tools import rename_species

import matplotlib as mpl

root = '/Users/xianggao/GPS data/prj_DNS'

mech = 'GRI'
t = 20
npz = os.path.join(root, mech,'detailed','raw','[syngas] + [air]','DNS','phi1.0_1.0atm_500.0K',mech+str(t)+'tj2D','raw.npz')
raw0 = load_raw(npz)

mech = 'SKE'
npz = os.path.join(root, mech,'detailed','raw','[syngas] + [air]','DNS','phi1.0_1.0atm_500.0K',mech+str(t)+'tj2D','raw.npz')
raw1 = load_raw(npz)

raw = dict()
for k in ['axis0','axis1','axis2']:
	raw[k] = raw0[k]
for k in ['heat_release_rate','temperature','mixture fraction']:
	raw[k] = raw0[k] - raw1[k]

mech = 'GRI-SKE'
fld_npz = os.path.join(root, mech,'detailed','raw','[syngas] + [air]','DNS','phi1.0_1.0atm_500.0K',mech+str(t)+'tj2D')
if not os.path.exists(fld_npz):
	os.makedirs(fld_npz)
npz = os.path.join(fld_npz,'raw.npz')
save_raw_npz(raw, npz)
