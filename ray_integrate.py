from radtrans_integrate import radtrans_integrate
import numpy as np
import scipy.integrate
import read_geodebug_file

# use radtrans_integrate LSODA to do rad trans along ray
# JAD 10/19/2016
def get_jar_geodebug(d):
    jj=np.zeros((len(d.ji),4))
    aa=np.zeros((len(d.ji),4))
    rr=np.zeros((len(d.ji),3))
    jj[:,0]=d.ji; jj[:,1]=d.jq; jj[:,2]=d.ju; jj[:,3]=d.jv
    aa[:,0]=d.ki; aa[:,1]=d.kq; aa[:,2]=d.ku; aa[:,3]=d.kv
    rr[:,0]=d.rhoq; rr[:,1]=d.rhou; rr[:,2]=d.rhov
    xx=d.lam
    return xx,jj,aa,rr

def run_geodebug(atol=1e-8,rtol=1e-6,sphstokes=-1,file='/Users/jdexter/code/grtrans/public/grtrans/geodebug.out'):
    d = read_geodebug_file.geodebug()
    d.file = file
    d.read()
    xx,jj,aa,rr = get_jar_geodebug(d)
    i = run(xx[::-1],jj,aa,rr,atol=atol,rtol=rtol,sphstokes=sphstokes)
    return i

def run(x,jarr,aarr,rhoarr,sphstokes=-1,atol=1e-8,rtol=1e-6):
    if sphstokes==-1:
        method=0
    else:
        method=3
    radtrans_integrate.init_radtrans_integrate_data(method,4,len(x),len(x),10.,0.1,atol,rtol,1e-2,100000)
    Karr = (np.append(aarr,rhoarr,axis=1))
    tau = np.append(0.,scipy.integrate.cumtrapz(Karr[:,0],x))
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    i = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    return i
