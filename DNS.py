# -*- coding: utf-8 -*-

from src.tools.xls2npz import *
from src.ct.def_ct_tools import *
import matplotlib.pyplot as plt
from src.tools.def_DNS import *
from src.core.def_GPSA import *
from src.post.def_painter import plot_GPedge, plot_GPedge_mf
from src.post.def_plt_tools import rename_species

import matplotlib as mpl





def convert_DNS_data(overwrite=False):

	for mech in ['GRI']:#['SKEGRItherm']:#,'SKE']:
		for case in ['0tj','10tj','20tj','30tj','40tj']:#,'30tj']:# 'ign'

			name = mech+case+'2D'

			#prj_dir = os.path.join('prj','prj_DNS','GRI')
			prj_dir = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)

			cti = os.path.join(prj_dir,'detailed','mech','chem.cti')
			soln = ct.Solution(cti)

			#sp = 'N2'; key = ('mole_fraction',soln.species_names.index(sp)); s_key = key[0]+'_'+sp
			key = 'active reactions'; s_key = key


			f_dat = os.path.join(prj_dir,name+'.dat')
			dir_raw = os.path.join(prj_dir,'detailed','raw','[syngas] + [air]',
				'DNS','phi1.0_1.0atm_500.0K',name)

			if not os.path.exists(dir_raw):
				os.makedirs(dir_raw)

			f_npz = os.path.join(dir_raw,'raw.npz')

			if not os.path.exists(f_npz) or overwrite:
				raw = DNS2raw(f_dat, soln)
				raw = save_raw_npz(raw, f_npz)
			else:
				print 'already exists: '
				print f_npz


def draw_DNS_contour():
	Z_st = 0.422

	plt.rc('font', **{'family':'Times New Roman'})

	cases = [str(i)+'tj' for i in [20,40]]
	#cases = [str(i)+'tj' for i in [20,20]]
	n_case = len(cases)

	key_x = 'mixture fraction'; xlim = (0,1); xlabel = 'Z'
	key_y = 'chi'; ylim = (-9,10); ylabel = 'chi'
	#key_y = 'temperature'; ylim = (400,2000); ylabel = 'T [K]'
	#key_y = 'heat_release_rate';  ylabel = 'Qdot (J/m3-s)'; ylim=(0,3e9)	# J/m3-s
 
	#key_x = 'axis0'; xlabel = 'x [mm]'; xlim=None
	#key_y = 'axis1'; ylabel = 'y [mm]'; ylim=None

	#sp = 'CO2'; ylim = (0,0.2); ylabel = sp
	

	#s_coord = 'TZ_'; 
	fig_h = 5
	s_coord = '_'.join([xlabel.split(' ')[0], ylabel.split(' ')[0]])


	for mech in ['GRI']:# ['GRIred','SKEGRItherm']:
		prj_dir = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)
		cti = os.path.join(prj_dir,'detailed','mech','chem.cti')
		soln = ct.Solution(cti)

		#key_y = ('mole_fraction',soln.species_names.index(sp));

		clim = None
		mutiplier = 1.0
		s_key = None

		#key = 'heat_release_rate'; s_key = 'Qdot (J/m3-s)'; clim=(-2.5e9,2.5e9)#(0,3e9)	# J/m3-s
		#key = 'cp'; s_key = 'cp'; clim=None#(-2.5e9,2.5e9)#(0,3e9)	# J/m3-s
		#key = 'pressure'; s_key = 'p'
		#key = 'temperature'; s_key = 'T [K]'; clim=(400,2000)
		#key = 'chi'; s_key = 'chi'; clim=(0, 10000)
		#key = 'mixture fraction'; s_key = 'MixFrac'
		

		#sp = 'HO2'; clim = (0,350); mutiplier = 1e6; s_key = sp+' (ppm)'
		#sp = 'OH'; clim = (0,0.01);s_key = sp
		#sp = 'OH'; clim = (-24,24)
		#sp = 'H'; clim = (0,0.03);
		#sp = 'O'; clim = (0,0.015);s_key = sp
		sp = 'H2O'; clim = (0,0.05);s_key = sp
		#sp = 'H2'; clim = (0,0.05);s_key = sp
		#sp = 'CO2'; clim = (0,0.2);s_key = sp
		key = ('mole_fraction',soln.species_names.index(sp));
		s_key = sp+';'+key[0]
		GPSA_data_list = None; Z_cut = None; 
		#Z_contour = (
		#			[0.05, Z_st],
		#			['Z=0.05','Z=Z'+r'$_{st}$']
		#			)
		Z_contour = None
		cmap = 'jet'
		#cmap = 'bwr'

		prj_json = os.path.join(prj_dir,'project.json')
		project = json.load(open(prj_json,'r'))




		fig, axs = plt.subplots(1,n_case,figsize=(5.5*n_case, fig_h))

		vmin = float('inf')
		vmax = -float('inf')
		for i_case in range(n_case):
			case = cases[i_case]

			name = mech+case+'2D'

			dir_raw = os.path.join(prj_dir,'detailed','raw','[syngas] + [air]',
				'DNS','phi1.0_1.0atm_500.0K',name)
			f_npz = os.path.join(dir_raw,'raw.npz')

			raw = load_raw(f_npz)

			""" -------
			key = 'GPSA'; 
			#method = 'D_GP'; clim = (0,1);#Z_cut = 0.05; 
			method = 'R_GP'; clim = (-2000,2000); cmap='bwr'; Z_contour=None
			
			traced = 'C'
			#GP_alias = 'GP-H2-HO2'#,
			#GP_alias = 'GP-H2-highT'
			GP_alias = 'CO --> CO2'
			s_key = method+';'+GP_alias
			GP_dir = None
			for GP_name in project['GP_'+traced].keys():
				if project['GP_'+traced][GP_name]['alias'] == GP_alias:
					GP_dir = project['GP_'+traced][GP_name]
					#print GP_dir
					break
			if GP_dir is None:
				print 'could not find GP whose alias is '+str(GP_alias)
				sys.exit()
			GPSA_data_list = load_GPSA(dir_raw, GP_dir, method)
			
			# ----- """ 


			print 'name = '+name
			print 'min(T) = '+str(min(raw['temperature']))
			ax = axs[i_case]
			im = show_2D(raw, key, GPSA_data_list=GPSA_data_list,
				key_x=key_x, key_y=key_y, xlim=xlim, ylim=ylim, cmap=cmap,
				clim=clim, ax=ax, mutiplier=mutiplier, Z_cut=Z_cut, Z_contour=Z_contour)

			ax.set_title(name+'\n')
			#ax.get_xaxis().set_visible(False)
			ax.set_xlabel(xlabel)
			if ylim is not None and i_case>0:
				ax.get_yaxis().set_visible(False)
			
			if i_case == 0:
				ax.set_ylabel(ylabel)
			else:
				ax.get_yaxis().set_visible(False)

			#xticks_old = ax.get_xticks()
			#print 'xticks_old = ', xticks_old
			#sys.exit()


			cax = None
			if clim is None:
				cax, kw = mpl.colorbar.make_axes([axs[i_case]], location='right')	
			elif i_case+1 == n_case:
				cax, kw = mpl.colorbar.make_axes(list(axs), location='right')	

			

			if cax is not None:
				cb = plt.colorbar(im, cax=cax, **kw)
				cb.set_label('\n'+str(s_key)+', '+mech)




		
		fld_plot = os.path.join(prj_dir,'plot','DNS raw')

		if not os.path.exists(fld_plot):
			os.makedirs(fld_plot)

		fig_name = s_coord +'_' + s_key.split('(')[0].strip()
		print fig_name
		
		#plt.subplots_adjust(top=0.97, bottom=0.5)

		path = os.path.join(fld_plot, fig_name+'.pdf')
		plt.savefig(path)
		print 'path = ',path
		plt.close()









def draw_DNS_edge_pnt(task='edge'):

	plt.rc('font', **{'family':'Times New Roman'})
	cases = [str(i)+'tj' for i in [20,40]]
	n_case = len(cases)

	target_x = 4.0
	target_y = 8.0

	mechs = ['GRIred','SKEGRItherm']

	if task == 'mf':
		f, axs = plt.subplots(1,2)


	for i_mech in range(len(mechs)):
		mech = mechs[i_mech]
		color_mech = ['b','r'][i_mech]

		prj_dir = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)
		cti = os.path.join(prj_dir,'detailed','mech','chem.cti')
		soln = ct.Solution(cti)

		prj_json = os.path.join(prj_dir,'project.json')
		project = json.load(open(prj_json,'r'))

		for i_case in range(n_case):
			case = cases[i_case]
			name = mech+case+'2D'

			dir_raw = os.path.join(prj_dir,'detailed','raw','[syngas] + [air]',
				'DNS','phi1.0_1.0atm_500.0K',name)
			f_npz = os.path.join(dir_raw,'raw.npz')
			raw = load_raw(f_npz)

			xx = raw['axis0']
			yy = raw['axis1']

			minx = min(xx)
			miny = min(yy)

			xx = [x-minx for x in xx]
			yy = [y-miny for y in yy]

			i_plot = None
			for i in range(len(xx)):
				if xx[i]>=target_x/1e3 and yy[i]>=target_y/1e3:
					i_plot = i
					T_plot = raw['temperature'][i_plot]
					print 'found i_plot = '+str(i_plot)+', x = '+str(xx[i_plot])+\
						', y = '+str(yy[i_plot])+', T = '+str(T_plot)
					break

			if i_plot is None:
				print 'cound not find i_plot, max(xx) = '+str(max(xx))+', max(yy) = '+str(max(yy))
				sys.exit()


			if task == 'edge':
				# ================================================

				key = 'GPSA'; 
				#method = 'D_GP'; clim = (0,1);Z_cut = 0.05; 
				method = 'R_ij';# clim = (-4000,4000); cmap='bwr'; Z_contour=None
				
				traced = 'H'
				#GP_alias = 'GP-H2-HO2'#,
				GP_alias = 'GP-H2-highT'
				#GP_alias = 'CO --> CO2'
				s_key = method+';'+GP_alias
				GP_dir = None
				for GP_name in project['GP_'+traced].keys():
					if project['GP_'+traced][GP_name]['alias'] == GP_alias:
						GP_dir = project['GP_'+traced][GP_name]
						#print GP_dir
						break
				if GP_dir is None:
					print 'could not find GP whose alias is '+str(GP_alias)
					sys.exit()

				GPSA_data_list = load_GPSA(dir_raw, GP_dir, method)
				
				fld_plot = os.path.join(prj_dir,'plot','DNS raw')
				fig_name = 'edge_'+GP_alias+'_'+str(case)+str((target_x,target_y))
				path_save = os.path.join(fld_plot, fig_name+'.pdf')

				opt = dict()
				opt['fig_w'] = [6]
				opt['fig_h'] = [9]
				opt['xscale'] = 'linear'
				opt['xlim'] = [0,4000]
				opt['n_rxn'] = [2]
				opt['method'] = method

				plot_GPedge(soln, GPSA_data_list, opt, raw, path_save, project['rename'], i_plot, 
					fig_name+', '+mech+', T = '+str(T_plot))
				#plot_GPedge_mf(soln, GP_dir, opt, raw, path_save, project['rename'], i_plot, fig_name+', '+mech)

				print path_save
				plt.close()
				#sys.exit()

			else:				
				# ================================================

				
				
				sp = soln.species_names
				mf = raw['mole_fraction']
				for i_sp in range(len(sp)):
					x = mf[i_plot, i_sp]
					y = 1.0*i_sp + 1.0*i_mech/len(mechs)
					axs[i_case].plot([1e-9,x],[y,y],color=color_mech)
					axs[i_case].text(x, y, sp[i_sp]+','+mech,color=color_mech)
					axs[i_case].set_xlim(0,0.1)

	if task == 'mf':
		fld_plot = os.path.join(prj_dir,'plot','DNS raw')
		fig_name = 'mf_'+str((target_x,target_y))
		path_save = os.path.join(fld_plot, fig_name+'.pdf')
		plt.savefig(path_save)
















def draw_DNS_stat():
	Z_st = 0.422

	plt.rc('font', **{'family':'Times New Roman'})

	#cases = [str(i)+'tj' for i in [20,40]]
	#tt = [0,10,20,30,40]
	tt = [20,40]
	cases = [str(i)+'tj' for i in tt]
	n_case = len(cases)

	n_x = 100

	key_x = 'mixture fraction'; xlim = (0,1); xlabel = 'Z'
	#key_x = 'chi'; xlim = (0,2e4); xlabel = 'chi'
	#key_y = 'temperature'; ylim = (400,2000); ylabel = 'T (K)'
	#key_y = 'heat_release_rate';  ylabel = 'Qdot (J/m3-s)'; ylim=(0,3e9)	# J/m3-s
 
	#key_x = 'axis0'; xlabel = 'x [mm]'; xlim=None
	#key_y = 'axis1'; ylabel = 'y [mm]'; ylim=None

	#sp = 'CO2'; ylim = (0,0.2); ylabel = sp
	

	#s_coord = 'TZ_'; 
	fig_h = 5
	#s_coord = '_'.join([xlabel.split(' ')[0], ylabel.split(' ')[0]])


	fig, axs = plt.subplots(1,n_case,figsize=(5.5*n_case, fig_h),sharey='row', sharex=True)

	mechs = ['GRI','SKE']#,'SKEGRItherm']
	#mechs = ['GRI','SKE']
	#mech_labels = ['GRI-Mech 3.0','11-Species Model ']
	mech_labels = mechs

	y_st = dict()

	for i_mech in range(len(mechs)):
		mech = mechs[i_mech]
		y_st[mech] = []

		color = ['r','b','k'][i_mech]
		ls = ['-','--','-.'][i_mech]

		prj_dir = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)
		cti = os.path.join(prj_dir,'detailed','mech','chem.cti')
		soln = ct.Solution(cti)

		#key_y = ('mole_fraction',soln.species_names.index(sp));

		clim = None
		mutiplier = 1.0
		s_key = None

		#key = 'heat_release_rate'; s_key = 'Qdot (J/m3-s)'; clim=(0,3e9)	# J/m3-s
		#key = 'pressure'; s_key = 'p'
		#key = 'temperature'; s_key = 'T [K]'; clim=(400,2000)
		#key = 'chi'; s_key = 'chi'; clim=(0, 10000)
		#key = 'mixture fraction'; s_key = 'MixFrac'
		

		sp = 'HO2'; clim = (0,350); mutiplier = 1e6; s_key = sp+' (ppm)'
		#sp = 'OH'; clim = (0,0.01);s_key = sp
		#sp = 'H'; clim = (0,0.03);s_key = sp
		#sp = 'O'; clim = (0,0.015);s_key = sp
		#sp = 'H2O'; clim = (0,0.05);s_key = sp
		#sp = 'H2'; clim = (0,0.05);s_key = sp
		#sp = 'CO2'; clim = (0,0.2);s_key = sp
		key = ('mole_fraction',soln.species_names.index(sp));
		GPSA_data_list = None; Z_cut = None; 
		#Z_contour = (
		#			[0.05, Z_st],
		#			['Z=0.05','Z=Z'+r'$_{st}$']
		#			)
		Z_contour = None
		cmap = 'jet'

		prj_json = os.path.join(prj_dir,'project.json')
		project = json.load(open(prj_json,'r'))

		vmin = float('inf')
		vmax = -float('inf')
		for i_case in range(n_case):
			case = cases[i_case]

			name = mech+case+'2D'

			dir_raw = os.path.join(prj_dir,'detailed','raw','[syngas] + [air]',
				'DNS','phi1.0_1.0atm_500.0K',name)
			f_npz = os.path.join(dir_raw,'raw.npz')

			raw = load_raw(f_npz)

			""" -------
			key = 'GPSA'; 
			#method = 'D_GP'; clim = (0,1); Z_contour=None; ylabel = r'$D_{GP}$'+' [-]'
			
			method = 'R_GP'; clim = (-4000,4000); cmap='bwr'; Z_contour=None; ylabel = r'$R_{GP}$'+' [mole/m'+r'$^3$'+'-s]'
			
			
			#GP_alias = 'GP-H2-HO2'; traced = 'H'
			GP_alias = 'GP-H2-highT'; traced = 'H'
			#GP_alias = 'CO --> CO2'; traced = 'C'
			s_key = method+';'+GP_alias
			
			GP_dir = None
			for GP_name in project['GP_'+traced].keys():
				if project['GP_'+traced][GP_name]['alias'] == GP_alias:
					GP_dir = project['GP_'+traced][GP_name]
					#print GP_dir
					break
			if GP_dir is None:
				print 'could not find GP whose alias is '+str(GP_alias)
				sys.exit()
			GPSA_data_list = load_GPSA(dir_raw, GP_dir, method)
			y = GPSA_data_list[:-1]; s_key
			#"""

			"""
			sp = 'O2'
			i_sp = soln.species_names.index(sp)
			y = raw['mole_fraction'][:,i_sp]; s_key = sp; ylabel = rename_species(sp, project['rename'])
			#"""

			#y = raw['heat_release_rate']; s_key = 'Qdot'; ylabel = 'Qdot'
			y = raw['cp']; s_key = 'cp'; ylabel = 'cp'
			#y = raw['chi']; s_key = 'chi'; ylabel = r'$\Chi$'
			#y = raw['temperature']; s_key = 'T'; ylabel = 'T [K]'

			x, xx, ii = find_coord(raw, key_x, max_len=n_x, lim=xlim)


			ygroup = [None]*len(xx)
			for j in range(len(y)):
				i = ii[j]
				if ygroup[i] is None:
					ygroup[i] = []
				ygroup[i].append(y[j])


			mean_y = []
			med_y = []
			max_y = []
			std_y = []
			x_show = []
			for i in range(n_x):
				if ygroup[i] is not None:
					mean_y.append(np.mean(ygroup[i]))
					med_y.append(np.median(ygroup[i]))
					max_y.append(max(ygroup[i]))
					std_y.append(np.std(ygroup[i]))

					x_show.append(xx[i])

			if xlabel == 'Z':
				diff = [abs(x-Z_st) for x in x_show]
				i_st = diff.index(min(diff))
				print 'find Zst at ',x_show[i_st]
				y_st[mech].append(mean_y[i_st])



			#ax = axs[1,i_case]
			#ax.plot(x_show,std_y,color=color, label=mech, linestyle=ls)#, marker='o', fillstyle='none')
			ax = axs[i_case]
			ax.plot(x_show,mean_y,color=color, label=mech, linestyle=ls)#, marker='o', fillstyle='none')
			#ax.plot(xx,[mean_y[j] + std_y[j] for j in range(n_x)],color=color,linestyle='--')#,linestyle='None',marker='.')
			#ax.plot(xx,[mean_y[j] - std_y[j] for j in range(n_x)],color=color,linestyle='--')#,linestyle='None',marker='.')

			if i_mech == 0:
				ax.set_title(ylabel+','+case+'\n\n')
				ax.set_xlabel(xlabel)
				#axs[1,i_case].set_xlabel(xlabel)
				#ax.set_ylim(clim)
				ax.set_xlim(xlim)
				ax.axhline(y=0, color='k', linestyle=':')


			if i_case == 0:
				ax.set_ylabel('mean of '+ylabel)
				#axs[1,0].set_ylabel('std of '+ylabel)
			#else:
			#	ax.get_yaxis().set_visible(False)



	fld_plot = os.path.join(prj_dir,'plot','DNS stat')
	if not os.path.exists(fld_plot):
		os.makedirs(fld_plot)
	fig_name = 'stat_'+xlabel+'_'+s_key
	axs[0].legend(frameon=False, loc='lower center')
	plt.tight_layout()
	path = os.path.join(fld_plot, fig_name+'.pdf')
	plt.savefig(path)
	#print fld_plot
	print 'path = ',path
	plt.close()


	if xlabel == 'Z':
		f, ax = plt.subplots(figsize=(6,3))


		for i_mech in range(len(mechs)):
			mech = mechs[i_mech]
			color = ['r','b','k'][i_mech]
			ls = ['-','--','-.'][i_mech]

			ax.plot(tt, y_st[mech],color=color,linestyle=ls,label=mech_labels[i_mech])
			print mech
			print tt
			print y_st[mech]

		ax.legend(frameon=False, loc='lower center')
		ax.set_xlabel('Time ['+r'$t_j$'+']')
		ax.set_ylabel('HRR(Z=Z'+r'$_{st}$'+') [J/m'+r'$^3$'+'-s]')

		ticks_old = ax.get_yticks()

		ticks = []
		for i in range(len(ticks_old)):
			if int(i/2) == 0.5*i:
				ticks.append(ticks_old[i])

		ticklabels = []
		s = r'$\times 10^9$'
		for tick in ticks:
			ticklabels.append(str(tick/1e9)+s)

		ax.set_ylim(0.8e9, 2.6e9)

		ax.set_yticks(ticks)
		ax.set_yticklabels(ticklabels)
		plt.tight_layout()
		plt.savefig(os.path.join(fld_plot, fig_name+' vs t lab'+'.pdf'))
		plt.close()















def plot_scatter(x0=None, y0=None, sp=None):
	plt.rc('font', **{'family':'Times New Roman'})

	ZZ=[0.31,0.422,0.54]
	dZ=0.01


	mech = 'GRI'

	prj_dir = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)
	cti = os.path.join(prj_dir,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)

	prj_json = os.path.join(prj_dir,'project.json')
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
		dir_raw = os.path.join(prj_dir,'detailed','raw','[syngas] + [air]',
				'DNS','phi1.0_1.0atm_500.0K',name)
		f_npz = os.path.join(dir_raw,'raw.npz')

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
		
		GP_dir = None
		for GP_name in project['GP_'+traced].keys():
			if project['GP_'+traced][GP_name]['alias'] == GP_alias:
				GP_dir = project['GP_'+traced][GP_name]
				#print GP_dir
				break
		if GP_dir is None:
			print 'could not find GP whose alias is '+str(GP_alias)
			sys.exit()
		GPSA_data_list = load_GPSA(dir_raw, GP_dir, method)
		y_all = GPSA_data_list
		#"""

		#"""
		#sp = 'OH'
		i_sp = soln.species_names.index(sp)
		y_all = [float(yi) for yi in raw['mole_fraction'][:,i_sp]]; s_key_y = sp; ylabel = rename_species(sp, project['rename'])
		#"""
		
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





	fld_plot = os.path.join(prj_dir,'plot','DNS scatter')
	if not os.path.exists(fld_plot):
		os.makedirs(fld_plot)
	fig_name = 'scatter_'+s_key_x+'_'+s_key_y
	axs[0].legend(frameon=False, loc='upper left')
	plt.tight_layout()
	path = os.path.join(fld_plot, fig_name+'.pdf')
	plt.savefig(path)
	print path





def draw_0tj(sp):
	plt.rc('font', **{'family':'Times New Roman'})

	ZZ=[0.2,0.422,0.6]
	dZ=0.01


	mech = 'GRI'

	prj_dir = os.path.join('/Users','xianggao','GPS data','prj_DNS',mech)
	cti = os.path.join(prj_dir,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)

	prj_json = os.path.join(prj_dir,'project.json')
	project = json.load(open(prj_json,'r'))


	fig_h = 5

	fig, ax = plt.subplots(1,1,figsize=(5.5, fig_h))#,sharey='row', sharex=True)

	name = mech+'0tj2D'
	print name
	dir_raw = os.path.join(prj_dir,'detailed','raw','[syngas] + [air]',
			'DNS','phi1.0_1.0atm_500.0K',name)
	f_npz = os.path.join(dir_raw,'raw.npz')

	raw = load_raw(f_npz)
	#x_all = raw['axis1']; xlabel = 'y [mm]'; s_key_x = 'y'
	x_all = raw['temperature']; xlabel = 'T [K]'; s_key_x = 'T'

	axis0 = raw['axis0']
	fixed_axis0 = np.median(raw['axis0'])

	x = []
	y = []; 
	#ylabel = 'T [K]'; s_key_y = 'T'

	"""
	iC = soln.element_index('C')
	iH = soln.element_index('H')
	ylabel = 'C/H mass ratio [-]'
	s_key_y = 'C-H_mass_ratio'
	"""

	#sp = 'OH'
	i_sp = soln.species_names.index(sp)
	y_all = [float(yi) for yi in raw['mole_fraction'][:,i_sp]]; s_key_y = sp; ylabel = rename_species(sp, project['rename'])
		



	for i in range(len(x_all)):
		if axis0[i] == fixed_axis0:

			x.append(x_all[i])
			y.append(y_all[i])
			
			#y.append(raw['mixture fraction'][i])
			#y.append(raw['temperature'][i])


			"""
			soln = raw2soln(soln, raw, i)
			emf_C = soln.elemental_mole_fraction(iC)
			emf_H = soln.elemental_mole_fraction(iH)

			if emf_H>0:
				y.append(emf_C/emf_H)
			else:
				y.append(float('nan'))
				"""




	#ax.plot(x,y,label=ylabel,marker='.')

	#plt.show()
	return x,y













if __name__ == '__main__':
	#convert_DNS_data()
	#draw_DNS_contour()
	#draw_DNS_edge_pnt()
	#draw_DNS_stat()

	sp = 'HO2'
	x,y = draw_0tj(sp)
	plot_scatter(x,y,sp)







