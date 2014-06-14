import sys
sys.path.insert(0,'/home/jbraden/Documents/python_scripts/')

import myfileutils as myf
import numpy as np
import matplotlib.pyplot as plt

infile="evolved_efunc.dat"
#infile="efunct_n512_L268.dat"

vval = 1.
kperp=0.35**0.5

period = np.pi*(1.+vval**2)**0.5/vval  # period of effective mass, not bg
vals=myf.read_timeblock(infile,[0,1,2,3,4])
#omega=kperp
omega=np.pi/period

tvals=[]
fvals=[]
pvals=[]
npartf=[]
npartp=[]
for i in range(len(vals)):
    tvals.append(vals[i][1])
    fvals.append(vals[i][2])
    pvals.append(vals[i][3])

tvals=np.array(tvals)
tvals=tvals.T
tvals = tvals / 2. / period  # normalize to breather period
fvals=np.array(fvals)
pvals=np.array(pvals)
npartf=fvals**2
npartp=pvals**2

npartx=0.5*(npartp + omega**2*npartf) / omega

npart=[]
for i in range(len(npartx)):
    npart.append(np.sum(npartx[i]))

plt.plot(tvals[0],np.log(npart), linewidth=2)
plt.ylabel(r'$\log(n_{\mathrm{eff}})$',fontsize=16)
plt.xlabel(r'$t/T_{\mathrm{breather}}$',fontsize=16)

#plt.xlim(0,8)
