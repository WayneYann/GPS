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
		method, _ = s_key.split(';')
		if method == 'D_GP':
			return r'$D_{GP}$'+' [-]', (0,1)
		elif method == 'R_GP':
			return r'$R_{GP}$'+' [mole/m'+r'$^3$'+'-s]', None
		

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
		return load_GPSA(fld_raw, GP_dict, method)
		"""
		


def find_data_raw(s_key, raw, soln=None, project=None, fld_raw=None):

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
		method, GP_alias = s_key.split(';')
		for GP_name in project['GP_'+traced].keys():
			if project['GP_'+traced][GP_name]['alias'] == GP_alias:
				GP_dict = project['GP_'+traced][GP_name]
				break
		
		return load_GPSA(fld_raw, GP_dict, method)[:-1]


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


def plot_OppDiff_scatter(s_key_y, s_key_x='Z', axs=None, color='r'):

	xlabel, xlim = find_label_lim(s_key_x)
	ylabel, ylim = find_label_lim(s_key_y)

	if axs is None:
		axs = get_axs()

	for i_case in range(n_case):
		ax = axs[i_case]
		case = cases[i_case]
		for d, ls in zip(dd, lsls):
			strain = float(U)/float(d)

			"""
			path = os.path.join(fld_prj, "OppDiff", "U"+str(U)+"_d"+d+"cm.csv")
			df = pd.read_csv(path, delimiter=',')
			x = find_data_chemkin(s_key_x, df)
			y = find_data_chemkin(s_key_y, df)
			"""

			name = "OppDiff_U"+U+"_d"+d+"cm"
			print name
			fld_raw = find_fld_raw(fld_prj, name)

			raw = load_raw(os.path.join(fld_raw,'raw.npz'))
			x = find_data_raw(s_key_x, raw, soln, project, fld_raw)
			y = find_data_raw(s_key_y, raw, soln, project, fld_raw)


			ax.plot(x,y, color=color, linestyle=ls, linewidth=2,
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

			

def plot_DNS_scatter(s_key_y, s_key_x, axs=None, 
	color='b', label="DNS", stat=False):

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
		X = find_data_raw(s_key_x, raw, soln, project, fld_raw)
		Y = find_data_raw(s_key_y, raw, soln, project, fld_raw)
		max_dX = (max(X)-min(X))/n_stat
		len_x = len(X)

		if i_case==0:
			ms = 10
		else:
			ms = 0.7

		if stat and i_case>0:
			sortedXY = [list(x) for x in zip(*sorted(zip(X, Y), key=lambda pair: pair[0]))]
			batch = len(X)/n_stat
			dx = (max(X) - min(X))/(n_stat-1)
			x = []
			y = []
			ey = []
			prev_i = 0
			while prev_i < len_x:
				batch_i = min(len_x - prev_i - 1, batch)
				while True:
					if batch_i < batch/10:
						break
					if sortedXY[0][prev_i+batch_i] - sortedXY[0][prev_i] > max_dX:
						batch_i /= 2
						#print "batch_i reduced to",batch_i
					else:
						break
					
				xx = sortedXY[0][prev_i: prev_i+batch_i]
				yy = sortedXY[1][prev_i: prev_i+batch_i]
				prev_i = prev_i+batch_i+1

				#xx = sortedXY[0][i*batch: (i+1)*batch]
				#yy = sortedXY[1][i*batch: (i+1)*batch]
				x.append(np.mean(xx))
				y.append(np.mean(yy))
				ey.append(np.std(yy))

			ax.plot(x, y, 
				linestyle='None', marker='.', color=color, markersize=10, label=label)
			
			# error bar
			alpha = 0.3
			for i in xrange(len(x)):
				# vertical
				ax.plot([x[i], x[i]], [y[i]+ey[i], y[i]-ey[i]], 
					color=color, alpha=alpha)

				# horizon
				for sign in [-1,1]:
					yi = y[i] + sign * ey[i]
					ax.plot([x[i]-max_dX/4, x[i]+max_dX/4], [yi, yi],
						color=color, alpha=alpha)


		else:
			ax.plot(X, Y, 
				linestyle='None', marker='.', color=color, markersize=ms, label=label)

		ax.set_title(case)
		ax.set_xlabel(xlabel)
		if xlim is not None:
			ax.set_xlim(xlim)
		if i_case==0:
			ax.set_ylabel(ylabel)
			if ylim is not None:
				ax.set_ylim(ylim)

	return axs



def plot_scatter_vs_OppDiff():

	s_key_x = 'Z'
	#s_key_x = 'T'
	stat = True

	for sp in ["CO2", "CO", "H2O"]:#['HO2', "H", "O", "OH",]:
		s_key_y = "sp;"+sp
	#for s_key_y in ['Qdot', 'T']:

		axs = plot_DNS_scatter(s_key_y, s_key_x, stat=stat)
		axs = plot_OppDiff_scatter(s_key_y, s_key_x, axs=axs)
		axs[0].legend(loc="best", frameon=False)

		fld_plot = os.path.join(fld_prj,'plot','scatter vs OppDiff')
		fig_name = s_key_y+'_'+s_key_x
		if stat:
			fig_name += '; stat'
		path = os.path.join(fld_plot, fig_name+'.png')
		makedirs(fld_plot)
		plt.savefig(path)
		plt.close()
		print path



def plot_scatter_GP():
	method = "R_GP"
	axs = None
	
	s_key_x = 'Z'
	#s_key_x = 'Qdot'
	#stat = True
	stat = False

	if stat:
		global dd
		dd = ["0.142"]
		global lsls
		lsls = ["-"]

	for GP_alias, color in zip(["GP-H2-HO2", "GP-H2-highT"], ['b','r']):
		s_key_y = method+";"+GP_alias
		axs = plot_DNS_scatter(s_key_y, s_key_x, 
			color=color, axs=axs, label=GP_alias, stat=stat)
		if stat:
			axs = plot_OppDiff_scatter(s_key_y, s_key_x, axs=axs, color=color)
	
	axs[0].legend(loc="best", frameon=False)
	fld_plot = os.path.join(fld_prj,'plot','scatter GP')
	fig_name = s_key_y+'_'+s_key_x
	if stat:
		fig_name += '; stat'
	path = os.path.join(fld_plot, fig_name+'.png')
	makedirs(fld_plot)
	plt.savefig(path)
	plt.close()
	print path




if __name__ == '__main__':
	# shape = 129 * 257 = 33153


	plt.rc('font', family='Times New Roman')

	n_stat = 15	# if plot mean, how many x levels?
	traced = 'H'
	mech = 'GRI'
	tt = [0,20,40]
	#ZZ=[0.31,0.422,0.54]	# list of mixture fraction to be plotted
	ZZ=[0.15,0.422,0.6]	# list of mixture fraction to be plotted
	dZ=0.005

	U = "69" 
	dd = ["0.1", "0.142", "0.9"]
	lsls = ["--",":","-"]

	cases = [str(i)+'tj' for i in tt]
	n_case = len(cases)
	fig_h = 5
	sub_w = 5.5
	fld_prj = find_fld_prj(mech)
	cti = os.path.join(fld_prj,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)
	prj_json = os.path.join(fld_prj,'project.json')
	project = json.load(open(prj_json,'r'))

	#plot_scatter_GP()
	plot_scatter_vs_OppDiff()



