# -*- coding: utf-8 -*-

from src.tools.xls2npz import *
from src.ct.def_ct_tools import *
import matplotlib.pyplot as plt
from src.tools.def_DNS import *
from src.core.def_GPSA import *
from src.post.def_painter import plot_GPedge, plot_GPedge_mf
from src.post.def_plt_tools import rename_species

import matplotlib as mpl
import pandas as pd


def find_label_lim(s_key):
	basics = {
		'Z'		: ('Z', (0,1)),
		'T'		: ('T [K]', (400, 2000)),
		'Qdot'	: ('HRR [J/m'+r'$^3$'+'-s]', None),
		}
	if s_key in basics:
		return basics[s_key]

	if 'sp' in s_key:
		_, sp = s_key.split(';')
		return sp, None

	if '_GP' in s_key:
		_, GP_alias = s_key.split(';')
		if method == 'D_GP':
			return r'$D_{GP}$'+' [-]'
		elif method == 'R_GP':
			return r'$R_{GP}$'+' [mole/m'+r'$^3$'+'-s]'
		


def find_data_chemkin(s_key, df, project=None):

	basics = {
		'Z'		: 'Mixture_fraction ()',
		'T'		: 'Temperature (K)',
		'Qdot'	: 'Net_heat_production_from_gas-phase_reactions (erg/cm3-sec)',
		}
	if s_key in basics:
		data = df[basics[s_key]]
		if s_key == 'Qdot':
			data = [x/10 for x in data]
		return data

	if 'sp' in s_key:
		_, sp = s_key.split(';')
		return df["Mole_fraction_"+sp+" ()"]

	if '_GP' in s_key:
		raise NotImplementedError
		"""
		method, GP_alias = s_key.split(';')
		for GP_name in project['GP_'+traced].keys():
			if project['GP_'+traced][GP_name]['alias'] == GP_alias:
				GP_dict = project['GP_'+traced][GP_name]
				break
				"""
		
		return load_GPSA(fld_raw, GP_dict, method)



def find_data_npz(s_key, raw, soln=None, project=None):

	basics = {
		'Z'		: 'mixture fraction',
		'T'		: 'temperature',
		'Qdot'	: 'heat_release_rate',
		}
	if s_key in basics:
		return raw[basics[s_key]]

	if 'sp' in s_key:
		_, sp = s_key.split(';')
		i_sp = soln.species_names.index(sp)
		return [float(yi) for yi in raw['mole_fraction'][:,i_sp]]

	if '_GP' in s_key:
		raise NotImplementedError
		"""
		method, GP_alias = s_key.split(';')
		for GP_name in project['GP_'+traced].keys():
			if project['GP_'+traced][GP_name]['alias'] == GP_alias:
				GP_dict = project['GP_'+traced][GP_name]
				break
				"""
		
		return load_GPSA(fld_raw, GP_dict, method)

"""
if method == 'D_GP':
			clim = (0,1)
			label = r'$D_{GP}$'+' [-]'

		elif method == 'R_GP':
			clim = (-4000,4000)
			cmap='bwr'
			label = r'$R_{GP}$'+' [mole/m'+r'$^3$'+'-s]'

"""



def find_fld_prj(mech):
	return os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)


def find_fld_raw(fld_prj, case_name):
	return os.path.join(fld_prj,
			'detailed','raw','[syngas] + [air]',
			'DNS','phi1.0_1.0atm_500.0K',
			case_name)


def makedirs(fld):
	if not os.path.exists(fld):
		os.makedirs(fld)


def get_axs():
	fig, axs = plt.subplots(1,n_case,
		figsize=(sub_w * n_case, fig_h),sharey='row', sharex=True)
	return axs


def plot_OppDiff_scatter(s_key_y, s_key_x='Z', axs=None):
	print 'ploting OppDiff scatter'

	xlabel, xlim = find_label_lim(s_key_x)
	ylabel, ylim = find_label_lim(s_key_y)

	if axs is None:
		axs = get_axs()

	U = "69" 
	dd = ["0.1", "0.142", "0.9"]
	lsls = ["--",":","-"]

	for i_case in range(n_case):
		ax = axs[i_case]
		case = cases[i_case]
		for d, ls in zip(dd, lsls):
			strain = float(U)/float(d)

			path = os.path.join(fld_prj, "OppDiff", "U"+str(U)+"_d"+d+"cm.csv")
			df = pd.read_csv(path, delimiter=',')
			x = find_data_chemkin(s_key_x, df)
			y = find_data_chemkin(s_key_y, df)

			ax.plot(x,y, color='r', linestyle=ls, linewidth=2,
				label="OppDiff, "+r"$\alpha=$"+str(int(strain))+r"$s^{-1}$")

		ax.set_title(case)
		ax.set_xlabel(xlabel)

		if xlim is not None:
			ax.set_xlim(xlim)
		if i_case==0:
			ax.set_ylabel(ylabel)
			if ylim is not None:
				ax.set_ylim(ylim)

	return axs

			

def plot_DNS_scatter(s_key_y, s_key_x='Z', axs=None):
	print 'ploting DNS scatter'

	xlabel, xlim = find_label_lim(s_key_x)
	ylabel, ylim = find_label_lim(s_key_y)

	if axs is None:
		axs = get_axs()

	for i_case in range(n_case):
		ax = axs[i_case]
		case = cases[i_case]
		name = mech+case+'2D'
		print name
		fld_raw = find_fld_raw(fld_prj, name)

		raw = load_raw(os.path.join(fld_raw,'raw.npz'))
		x = find_data_npz(s_key_x, raw, soln)
		y = find_data_npz(s_key_y, raw, soln)

		if i_case==0:
			ms = 10
		else:
			ms = 0.7
		ax.plot(x,y, linestyle='None', marker='.', color='b', markersize=ms, label="DNS")
		ax.set_title(case)
		ax.set_xlabel(xlabel)

		if xlim is not None:
			ax.set_xlim(xlim)
		if i_case==0:
			ax.set_ylabel(ylabel)
			if ylim is not None:
				ax.set_ylim(ylim)

	return axs





def plot_scatter_givenZ(x0=None, y0=None, sp=None):
	#plt.rc('font', **{'family':'Times New Roman'})
	plt.rc('font', family='Times New Roman')

	mech = 'GRI'

	fld_prj = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)
	cti = os.path.join(fld_prj,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)

	prj_json = os.path.join(fld_prj,'project.json')
	project = json.load(open(prj_json,'r'))


	tt = [0,20,40]
	cases = [str(i)+'tj' for i in tt]
	n_case = len(cases)

	fig_h = 5

	fig, axs = plt.subplots(1,n_case,figsize=(5.5*n_case, fig_h),sharey='row', sharex=True)

	iC = soln.element_index('C')
	iH = soln.element_index('H')



	#react = ['H2','CO']
	#prod = ['H2O','CO2']
	#radical = ['H','O','OH']

	#ii_react = [soln.species_index(react[i]) for i in range(len(react))]
	#ii_prod = [soln.species_index(prod[i]) for i in range(len(prod))]




	for i_case in range(n_case):
		case = cases[i_case]
		name = mech+case+'2D'
		print name
		fld_raw = os.path.join(fld_prj,'detailed','raw','[syngas] + [air]',
				'DNS','phi1.0_1.0atm_500.0K',name)
		f_npz = os.path.join(fld_raw,'raw.npz')

		raw = load_raw(f_npz)
		#print raw.keys()
		Z = raw['mixture fraction']
		x_all = raw['temperature']; s_key_x = 'T'; xlabel = 'T [K]'
		#x_all = raw['chi']; s_key_x = 'chi'; xlabel = 'chi'
		
		"""
		sp = 'H2O'
		i_sp = soln.species_names.index(sp)
		x_all = [float(xi) for xi in raw['mole_fraction'][:,i_sp]]; s_key_x = sp; xlabel = rename_species(sp, project['rename'])
		"""		

		""" -------
		key = 'GPSA'; 
		method = 'D_GP'; clim = (0,1); Z_contour=None; ylabel = r'$D_{GP}$'+' [-]'
		#method = 'R_GP'; clim = (-4000,4000); cmap='bwr'; Z_contour=None; ylabel = r'$R_{GP}$'+' [mole/m'+r'$^3$'+'-s]'
		
		
		GP_alias = 'GP-H2-HO2'; traced = 'H'
		#GP_alias = 'GP-H2-highT'; traced = 'H'
		#GP_alias = 'CO --> CO2'; traced = 'C'
		s_key_y = method+';'+GP_alias
		
		GP_dict = None
		for GP_name in project['GP_'+traced].keys():
			if project['GP_'+traced][GP_name]['alias'] == GP_alias:
				GP_dict = project['GP_'+traced][GP_name]
				#print GP_dict
				break
		if GP_dict is None:
			print 'could not find GP whose alias is '+str(GP_alias)
			sys.exit()
		GPSA_data_list = load_GPSA(fld_raw, GP_dict, method)
		y_all = GPSA_data_list
		#"""

		#"""
		i_sp = soln.species_names.index(sp)
		y_all = [float(yi) for yi in raw['mole_fraction'][:,i_sp]]; s_key_y = 'sp_'+sp; ylabel = rename_species(sp, project['rename'])
		#"""

		#y_all = raw['mixture fraction']; s_key_y = 'Z'; ylabel = 'Z'
		
		
		#y_all = raw['heat_release_rate']; s_key_y = 'Qdot'; ylabel = 'HRR [J/m'+r'$^3$'+'-s]'
		#y_all = raw['chi']; s_key_y = 'chi'; ylabel = 'chi'

		#ylabel = 'C/H mole ratio [-]'; s_key_y = 'C-H_mole_ratio'

		for i_Z in range(len(ZZ)):
			color = ['b','k','r'][i_Z]
			Z0 = ZZ[i_Z]

			x = []
			y = []
			for i in range(len(raw['axis0'])):
				if abs(Z[i]-Z0)<dZ:
					x.append(x_all[i])
					y.append(y_all[i])

					"""
					soln = raw2soln(soln, raw, i)
					emf_C = soln.elemental_mole_fraction(iC)
					emf_H = soln.elemental_mole_fraction(iH)
					y.append(emf_C/emf_H)
					"""

			#print 'pnt selected = '+str(len(y))
			ax = axs[i_case]

			if x0 is not None and i_Z==0:
				ax.plot(x0,y0,color='y',label='1D OppDiff')

			ax.plot(x,y,linestyle='None',marker='.',color=color,label='Z = '+str(Z0))

			ax.set_title(case)
			ax.set_xlabel(xlabel)
			if i_case==0:
				ax.set_ylabel(ylabel)





	fld_plot = os.path.join(fld_prj,'plot','DNS scatter')
	if not os.path.exists(fld_plot):
		os.makedirs(fld_plot)
	fig_name = 'scatterZ_'+s_key_x+'_'+s_key_y
	axs[0].legend(frameon=False, loc='upper left')
	plt.tight_layout()
	path = os.path.join(fld_plot, fig_name+'.pdf')
	plt.savefig(path)
	print path
	#plt.show()















if __name__ == '__main__':
	plt.rc('font', family='Times New Roman')

	mech = 'GRI'
	tt = [0,20,40]
	#ZZ=[0.31,0.422,0.54]	# list of mixture fraction to be plotted
	ZZ=[0.15,0.422,0.6]	# list of mixture fraction to be plotted
	dZ=0.005

	cases = [str(i)+'tj' for i in tt]
	n_case = len(cases)
	fig_h = 5
	sub_w = 5.5
	fld_prj = find_fld_prj(mech)
	cti = os.path.join(fld_prj,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)
	prj_json = os.path.join(fld_prj,'project.json')
	project = json.load(open(prj_json,'r'))

	s_key_x = 'Z'

	#s_key_y = "Qdot"
	#s_key_y = 'T'
	#for sp in ["H", "O", "OH", 'HO2']:
	#	s_key_y = "sp;"+sp

	for s_key_y in ["sp;HO2"]:

		axs = plot_DNS_scatter(s_key_y, s_key_x)
		axs = plot_OppDiff_scatter(s_key_y, s_key_x, axs=axs)

		axs[0].legend(loc="best", frameon=False)
		fld_plot = os.path.join(fld_prj,'plot','OppDiff vs DNS')
		fig_name = s_key_y+'_'+s_key_x
		path = os.path.join(fld_plot, fig_name+'.png')
		makedirs(fld_plot)
		plt.savefig(path)
		plt.close()
		print path


	#convert_DNS_data()
	#draw_DNS_contour()
	#draw_DNS_edge_pnt()
	#draw_DNS_stat()

	"""
	for sp in ['H','O','OH']:#,'H']:#['CO2', 'CO', 'OH', 'H2O']:
		x,y = draw_0tj(sp)
		plot_scatter_givenZ(x,y,sp)

		fld_plot = os.path.join(fld_prj,'plot','DNS scatter')
		fig_name = 'scatter_'+s_key_x+'_'+s_key_y
		path = os.path.join(fld_plot, fig_name+'.png')
		makedirs(fld_plot)
		plt.savefig(path)
		print path
		"""


	#for sp in ['H2O2']:#'OH','O','H']:
	#	plot_scatter('sp_'+sp)
	#for k in ['T','Qdot']:
	#	plot_scatter(k)


	#	GP_alias = 'GP-H2-HO2'; traced = 'H'
	#GP_alias = 'GP-H2-highT'; traced = 'H'







