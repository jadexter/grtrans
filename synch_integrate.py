from radtrans_integrate import radtrans_integrate
from polsynchemis import polsynchemis
import numpy as np
import scipy.integrate

# calculate synchrotron emissivity for given coefficients
def synch_jarho(nu,n,B,T,theta):
    if ((np.isscalar(nu)==False) & (np.isscalar(n)==True)):
        n = n + np.zeros(len(nu))
        B = B + np.zeros(len(nu))
        T = T + np.zeros(len(nu))
        theta = theta + np.zeros(len(nu))
    e = polsynchemis.polsynchth(nu,n,B,T,theta)
    j = e[:,:4]; a = e[:,4:8]; rho = e[:,8:]
    return j,a,rho

def run(x,jarr,aarr,rhoarr,sphstokes=-1,atol=1e-8,rtol=1e-6,max_tau=10):
    if sphstokes==-1:
        method=0
    else:
        method=3
    radtrans_integrate.init_radtrans_integrate_data(method,4,len(x),len(x),max_tau,0.1,atol,rtol,1e-2,100000)
    Karr = (np.append(aarr,rhoarr,axis=1))
    tau = np.append(0.,scipy.integrate.cumtrapz(Karr[:,0],x))
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    i = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    return i
