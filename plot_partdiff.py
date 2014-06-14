import numpy as np

infile = "evolved_particles.dat"
nlat=128
k2=0.05

a=np.genfromtxt(infile,usecols=[0])
x=np.reshape(a,(-1,nlat))
a=np.genfromtxt(infile,usecols=[1])
t=np.reshape(a,(-1,nlat))
a=np.genfromtxt(infile,usecols=[2])
phi=np.reshape(a,(-1,nlat))
a=np.genfromtxt(infile,usecols=[3])
dphi=np.reshape(a,(-1,nlat))

npart=k2*phi**2+dphi**2

np_left=[]
np_right=[]
np_tot=[]
for i in range(len(npart)):
    np_tot.append(np.sum(npart[i]))
    np_left.append(np.sum(npart[i][0:nlat/2-1]))
    np_right.append(np.sum(npart[i][nlat/2:nlat]))
