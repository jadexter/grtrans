import geokerr_interface
import numpy as np
from scipy.optimize import minimize

g = geokerr_interface.geokerr()

def cartesian(u,mu,phi):
    x = 1./u*np.sqrt(1.-mu**2.)*np.cos(phi)
    y = 1./u*np.sqrt(1.-mu**2.)*np.sin(phi)
    z = 1./u*mu
    return x,y,z

def calc_orbit(ab,args):
    a=args[0]; mu0=args[1]; uf=args[2]; tpr=args[3]
    u,mu,t,lam,phi,tpm,tpr,i1,i2 = g.run_camera(n=[1,1,2],avals=[ab[0],ab[0],ab[1],ab[1]],a=a,mu0=mu0,standard=1,mufill=0,kext=0,next=0,uf=uf,uout=1e-6,tpr=tpr,offset=0)
    x,y,z = cartesian(u[0,-1],mu[0,-1],-phi[0,-1])
#    print 'calc_orbit: ',a,mu0,uf
#    print 'calc_orbit: ',u,mu,phi
    return np.array([x,y,z]),np.array([u[0,-1],mu[0,-1],phi[0,-1],tpm[0,-1],tpr[0,-1],g.sm,g.su])

def calc_sep(xobs,xem):
   return np.sum((xobs[0]-xem[0])**2.+(xobs[1]-xem[1])**2.+(xobs[2]-xem[2])**2.)

def iterate(ab,args):
    xobs,robs = calc_orbit(ab,args[1:])
    sep = calc_sep(xobs,args[0])
#    print 'iterate: ',ab,xobs,sep
    return sep

# take first guess from my attempted flat space alpha, beta
def run(ustar,mustar,phistar,abguess,a=0.99,mu0=np.cos(45./180.*np.pi),tpr=-1):
    xs,ys,zs = cartesian(ustar,mustar,phistar)
    if tpr==-1:
        tpr=-(np.sign(ys)-1)/2
    print 'xs: ',xs,ys,zs,tpr
    args=[np.array([xs,ys,zs]),a,mu0,ustar,tpr]
    res = minimize(iterate,abguess,args=args)
    return res
