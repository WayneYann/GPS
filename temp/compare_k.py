
import os, sys, json
import cantera as ct
import numpy as np

root = '/Users/xianggao/GPS data/prj_DNS'
mechs = ['GRIred','SKE']

n_H2 = 10
n_CO = 50
n_O2 = n_H2 * 0.5 + n_CO * 0.5
n_N2 = n_H2 * 4 + n_O2 * 3

X = 'H2:'+str(n_H2)+','+\
	'CO:'+str(n_CO)+','+\
	'O2:'+str(n_O2)+','+\
	'N2:'+str(n_N2)

print X


T = 1000
P = 1.0 * ct.one_atm


rxns_comp = dict()
for mech in mechs:
	cti = os.path.join(root,mech,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)
	soln.TPX = T,P,X


	#print dir(soln)
	print soln.species_names[0]
	print soln.molecular_weights[0]
	#print len(soln.net_production_rates)
	sys.exit()
	

	for id_rxn in range(soln.n_reactions):


		rxn = soln.reaction(id_rxn)
		three_body = 'M' in str(rxn)

		#print rxn
		d = rxn.products
		sp_mu = rxn.reactants
		for sp in sp_mu.keys():
			if sp in d.keys():
				d[sp] -= sp_mu[sp]
				three_body = True
			else:
				d[sp] = -sp_mu[sp]

			if abs(d[sp]) == 0:
				del d[sp]
				

		sorted_sp = sorted(d.keys())

		name = ''
		name_sign = ''
		sign = int(np.sign(d[sorted_sp[0]]))


		for sp in sorted_sp:
			mu = d[sp]
			if int(mu)==mu:
				mu = int(mu)

			name += str(mu)+sp+';'
			name_sign += str(sign * mu)+sp+';'

		if three_body:
			name = '(three body)' + name

		
		if name not in rxns_comp.keys():
			rxns_comp[name] = dict()
			rxns_comp[name]['name_sign'] = name_sign

		s_id = str(id_rxn)
		kf = soln.forward_rate_constants[id_rxn]
		kr = soln.reverse_rate_constants[id_rxn]

		eqn = str(rxn)


		if mech in rxns_comp[name].keys():
			#rxns_comp[name][mech+'_duplicated'] = info
			print 'duplicated!!!!'
			rxns_comp[name][mech]['id'] += ';'+s_id
			rxns_comp[name][mech]['eqn'] += ';'+eqn
			rxns_comp[name][mech]['kf'] += kf
			rxns_comp[name][mech]['kr'] += kr

		else:
			rxns_comp[name][mech] = dict()
			rxns_comp[name][mech]['id'] = s_id
			rxns_comp[name][mech]['eqn'] = eqn
			rxns_comp[name][mech]['kf'] = kf
			rxns_comp[name][mech]['kr'] = kr




header = ['rxn_name','sign_name']
kk = ['id','eqn','kf','kr']

for key in kk:
	for mech in mechs:
		header.append(key+'; '+mech)


print header
with open('rxn_compare.csv','w') as f:
	f.write(','.join(header)+'\n')

	for name in rxns_comp.keys():
		ss = [name, rxns_comp[name]['name_sign']]
		for key in kk:
			for mech in mechs:
				if mech in rxns_comp[name].keys():
					ss.append(str(rxns_comp[name][mech][key]))
				else:
					ss.append('-')
		
		f.write(','.join(ss)+'\n')
		















"""


R = 8.314 * 1e3


for mech in mechs:
	print '-'*5
	cti = os.path.join(root,mech,'detailed','mech','chem.cti')
	soln = ct.Solution(cti)
	soln.TPX = T,P,X


	id_rxn = 7

	k = soln.forward_rate_constants[id_rxn]


	break


"""