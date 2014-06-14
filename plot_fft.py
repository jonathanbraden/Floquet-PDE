import numpy as np

#infile = "eigenfunction_1period_n1024_v1.dat"
infile = "evolved_efunc.dat"
#period=3.4004353847414772  # for 2**0.5-1
#floquet=0.36253624190137057 #for k^2=0.2
#floquet=floquet/period
#floquet=0.26323842940233555/period # for 2**0.5-1
#floquet=0.26323842940233089/period
#xrnge=(-50,50)

period=4.4428829381583661
#floquet=0.79792182285373026/period
#floquet=1.3577621007934497/period
floquet=1.2699527643261310/period
xrnge=(-20,20)

#period=2.*np.pi
#floquet=1.3886878591809975/period
#xrnge=(-30,30)

nlat=512

a=np.genfromtxt(infile,usecols=[0])
x=np.reshape(a,(-1,nlat))
a=np.genfromtxt(infile,usecols=[1])
t=np.reshape(a,(-1,nlat))
t2=t.T
a=np.genfromtxt(infile,usecols=[2])
f=np.reshape(a,(-1,nlat))
f2=f.T

fn=np.exp(-floquet*t)*f
fn2=fn.T
# x,t,f are now positions, times, and function values as a function of x in first column, versus time

# Now perform the FFTs
ntime=len(f)
ft=np.fft.rfft(fn,(ntime-1),axis=0)
ft2=ft.T

# Normalize the spectrum to the RMS at that location
flucpow=[]
fnorm=[]
fmax=[]
for i in range(len(ft2)):
    fpow=np.sum(abs(ft2[i])**2)
    flucpow.append(fpow)
    fnorm.append(ft2[i]/np.sqrt(fpow))
    fmax.append(np.max(abs(fnorm[i])))

fnorm=np.array(fnorm)
fnormt=fnorm.T
import matplotlib.pyplot as plt

# Now make whatever plots I want
plt.plot(x[0],abs(fnormt[1]),label=r'$\omega = 2\pi T_{\mathrm{breather}}^{-1}$')
plt.plot(x[0],abs(fnormt[3]),label=r'$\omega = 4\pi T_{\mathrm{breather}}^{-1}$')
plt.plot(x[0],abs(fnormt[5]),label=r'$\omega = 6\pi T_{\mathrm{breather}}^{-1}$')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$\frac{|\delta\tilde{\phi}_{\omega}(x)|}{\sigma_x}$',fontsize=16)
plt.xlim(xrnge)
plt.ylim(0,1.1)
plt.legend(loc='center')

plt.show()

# Before plotting the angles, turn them into continuous variables
theta=np.zeros([3,len(ft2)])
theta[0,:]=np.angle(ft2[:,1])
theta[1,:]=np.angle(ft2[:,3])
theta[2,:]=np.angle(ft2[:,5])

# A better plan here would be to do some interpolations to determine the new x's
thetanew=[]
xnew=[]
for i in range(len(theta)):
    pos=np.where(np.abs(np.diff(theta[i])) > np.pi)[0]
    thetanew.append(np.insert(theta[i],pos+1,np.nan))
    xnew.append(np.insert(x[0],pos+1,np.nan))

plt.plot(xnew[0],thetanew[0],label=r'$\omega = 2\pi T_{\mathrm{breather}}^{-1}$')
plt.plot(xnew[1],thetanew[1],label=r'$\omega = 4\pi T_{\mathrm{breather}}^{-1}$')
plt.plot(xnew[2],thetanew[2],label=r'$\omega = 6\pi T_{\mathrm{breather}}^{-1}$')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$\mathrm{arg}(\delta\tilde{\phi}_{\omega})$',fontsize=16)
plt.ylim(-np.pi,np.pi)
plt.xlim(xrnge)
plt.legend(loc='center')

plt.show()

def continuous_angles(y):
    mid=len(y)/2-1
    wind=np.zeros(len(y))
    d=np.diff(y)
    winding_num=0
    for i in range(len(d)):
        if (d[i]>np.pi):
            winding_num-=1
        elif (d[i]<-np.pi):
            winding_num+=1
        wind[i+1]=winding_num

    wind=wind-wind[mid]
    return y+2.*np.pi*wind

thetanew=[]
xnew=[]
plt.plot(x[0],continuous_angles(theta[0]),label=r'$\omega = 2\pi T_{\mathrm{breather}}^{-1}$')
plt.plot(x[0],continuous_angles(theta[1]),label=r'$\omega = 4\pi T_{\mathrm{breather}}^{-1}$')
plt.plot(x[0],continuous_angles(theta[2]),label=r'$\omega = 6\pi T_{\mathrm{breather}}^{-1}$')
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$\mathrm{arg}(\delta\tilde{\phi}_{\omega})$',fontsize=16)
plt.xlim(xrnge)
plt.ylim(-5*np.pi,np.pi)
plt.legend(loc='lower left')

# Define a generator to transform my angular variables
import itertools
def wrap_angles(x,y):
    val1, val2 = itertools.tee(values)
# Get the first element and step val2 ahead one
    yield next(val2)
    for (px, py), (cx, cy) in itertools.izip(val1,val2):
        if (abs(cy-py) > np.pi):
            if (cy > py):
                yield (cx, py-2.*np.pi)
                yield (cx, None)
                yield (px, py+2.*np.pi)
            elif (cy < py):
                yield (cx,cy+2.*np.pi)
                yield (cx,None)
                yield (px,py-2.*np.pi)
        yield (cx,cy)
