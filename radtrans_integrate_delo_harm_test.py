import numpy as np

rhoarr = np.zeros((len(rhov),3))
jarr = np.zeros((len(ji),4))
aarr = np.zeros((len(ki),4))
rhoarr[:,0] = rhoq
rhoarr[:,1] = rhou
rhoarr[:,2] = rhov

jarr[:,0] = ji
jarr[:,1] = jq
jarr[:,2] = ju
jarr[:,3] = jv

aarr[:,0] = ki
aarr[:,2] = ku
aarr[:,1] = kq
aarr[:,3] = kv

Karr = np.append(aarr,rhoarr,axis=1)
tau = np.append(0.,scipy.integrate.cumtrapz(Karr[:,0],lam[::-1]))
radtrans_integrate.init_radtrans_integrate_data(1,4,len(lam),len(lam),10.,0.1,1e-8,1e-6,1e-2)
radtrans_integrate.integrate(lam,jarr[:,:],Karr[:,:],tau,4)
intens = radtrans_integrate.intensity.copy()
