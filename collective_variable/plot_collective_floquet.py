import matplotlib.pyplot as plt
import numpy as np

files=["floquet_r2.5_w40.dat","floquet_r3_w40.dat","floquet_r3.5_w40.dat"]
labels=[r'$r_{max}=2.5\sqrt{2}$',r'$r_{max}=3\sqrt{2}$',r'$r_{max}=3.5\sqrt{2}$']
markers=['b.','g^','rv']

for i in range(len(files)):
    a=np.genfromtxt(files[i],usecols=[0,1,2])
    plt.plot(a[:,1]**0.5*a[:,0],a[:,2],markers[i],label=labels[i],linewidth=1.5, markersize=3)

plt.xlabel(r'$k_\perp T_{walls}$',fontsize=26)
plt.ylabel(r'$\mu_{max}T_{walls}$',fontsize=30)
plt.xlim(0,20)
plt.ylim(0,18)
plt.legend()
plt.savefig('floquet_2walls_max.png')
plt.show()

for i in range(len(files)):
    a=np.genfromtxt(files[i],usecols=[0,1,3])
    plt.plot(a[:,1]**0.5*a[:,0],a[:,2],markers[i][0],label=labels[i],linewidth=1.5)

plt.xlabel(r'$k_\perp T_{walls}$',fontsize=20)
plt.ylabel(r'$\mu_{max,2} T_{walls}$',fontsize=20)
plt.xlim(0,20)
plt.legend()
plt.savefig('floquet_2walls_max2.pdf')
plt.show()
