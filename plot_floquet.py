import numpy as np

#infile="floquet_sg_n64_rat1e-5.dat"
#infile="floquet_smallv_n64_sg.dat"
infile="floquet_sinegordon_spectral_n64.dat"
numk=201

a=np.genfromtxt(infile,usecols=[0,1,2])
vvals=a[:,0]
vvals=np.reshape(vvals,(-1,numk))
vvals=vvals.T

kvals=a[:,1]
kvals=np.reshape(kvals,(-1,numk))
kvals=kvals.T

muvals=a[:,2]
muvals=np.reshape(muvals,(-1,numk))
muvals=muvals.T

# Redefine variables into appropriate coordinates
kvals=kvals*(1.+vvals**2)/vvals**2

import matplotlib.pyplot as plt
plt.contourf(1./vvals,kvals,2.*muvals,31,cmap=plt.cm.OrRd)
cb=plt.colorbar(orientation='vertical')
cb.ax.set_ylabel(r'$[Re(\mu)_{\mathrm{max}}]T_{\mathrm{breather}}$',fontsize=16)
plt.ylabel(r'$k_\perp^2(1+v^{-2})$',fontsize=16)
plt.xlabel(r'$v^{-1}$',fontsize=16)

plt.show()
