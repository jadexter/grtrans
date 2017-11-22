import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
plt.ion()
from radtrans_integrate import radtrans_integrate
rhoarr = np.zeros((300,3))
jarr = np.tile(np.array([1.,0.,0.,0.]),300).reshape(300,4)
a = np.array([10.,4.,1.,7.])
a2 = np.array([5.,3.,2.,1.])
a2arr = np.tile(a2,150).reshape(150,4)
aarr = np.tile(a,150).reshape(150,4)
aarr = np.append(a2arr/3.,3.*aarr,axis=0)
Karr = (np.append(aarr,rhoarr,axis=1))
x = np.cumsum(np.zeros(300)+1e-2)
tau = np.append(0.,scipy.integrate.cumtrapz(Karr[:,0],x))
radtrans_integrate.init_radtrans_integrate_data(1,4,300,300,10.,0.1,1e-8,1e-6,1e-2)
radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)

