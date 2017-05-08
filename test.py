import cantera as ct
import os
import matplotlib.pyplot as plt
from scipy.stats.mstats import gmean
import numpy as np

"""
#cti = os.path.join('prj','prj_multi-sL','detailed','mech','chem.cti')
soln = ct.Solution('gri30.xml')
rxn = soln.reaction_equations()

print sum(soln.delta_enthalpy * soln.net_rates_of_progress)
#"""

print np.arange(0,1,0.1)