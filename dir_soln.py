import cantera as ct

soln = ct.Solution('gri30.xml')
#d = '\n'.join([str(k) for k in dir(soln)])
#print d

T = 300
P = ct.one_atm
X = 'CH4:1'

soln.TPX = T,P,X

print soln.element_names
print soln.element_index('C')
#print soln.elemental_mole_fraction(2)
