import matplotlib.pyplot as plt
import numpy as np

AbE = [7.497000e+17, -1.41, 29580.07]

A = AbE[0]
b = AbE[1]
E = AbE[2] * 4.184	# J/mol

R = 8.314 # J/mol-K

TT = range(600,900,20)
kk = []
for T in TT:
	kk.append(A * (T**b) * np.exp(-E/(R*T)))

#plt.semilogy(TT, kk)
#plt.show()

print 'c7h15-2 => pc4h9 + c3h6'.upper()