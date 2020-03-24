from __future__ import print_function
import numpy as np
#import numpy.linalg
import scipy.integrate

# solve polarized RT equation analytically either using matricant (O-matrix) method from Degl'Innocenti or DELO method from Rees+
# JAD 8/12/2014
def opacity_matrix(a,p):
    return np.array([[a[0],a[1],a[2],a[3]],[a[1],a[0],p[2],-p[1]],[a[2],-p[2],a[0],p[0]],[a[3],p[1],-p[0],a[0]]])

def imatrix_4_test(a):

    a11 = a[0,0]; a12 = a[0,1]; a13 = a[0,2]; a14 = a[0,3]
    a21 = a[1,0]; a22 = a[1,1]; a23 = a[1,2]; a24 = a[1,3]
    a31 = a[2,0]; a32 = a[2,1]; a33 = a[2,2]; a34 = a[2,3]
    a41 = a[3,0]; a42 = a[3,1]; a43 = a[3,2]; a44 = a[3,3]

    b = np.zeros((4,4))
    
    b[0,0] = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42
    b[0,1] = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43
    b[0,2] = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42
    b[0,3] = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33
    b[1,0] = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43
    b[1,1] = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41
    b[1,2] = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43
    b[1,3] = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31
    b[2,0] = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41
    b[2,1] = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42
    b[2,2] = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41
    b[2,3] = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32
    b[3,0] = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42
    b[3,1] = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41
    b[3,2] = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42
    b[3,3] = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31
    detA = a[0,0]*b[0,0] + a[1,0]*b[0,1] +a[2,0]*b[0,2] + a[3,0]*b[0,3]
    b = b/detA
    return b,detA

# analytic inverse of a 4x4 matrix
def imatrix_4(a):
    a11 = a[0,0,:]; a12 = a[0,1,:]; a13 = a[0,2,:]; a14 = a[0,3,:]
    a21 = a[1,0,:]; a22 = a[1,1,:]; a23 = a[1,2,:]; a24 = a[1,3,:]
    a31 = a[2,0,:]; a32 = a[2,1,:]; a33 = a[2,2,:]; a34 = a[2,3,:]
    a41 = a[3,0,:]; a42 = a[3,1,:]; a43 = a[3,2,:]; a44 = a[3,3,:]
    a22a33a44 = a22*a33*a44; a23a34a42 = a23*a34*a42
    detA = a11*a22a33a44 + a11*a23a34a42 + a11*a24*a32*a43 + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42

    b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42
    b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43
    b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42
    b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33
    b21 = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43
    b22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41
    b23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43
    b24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31
    b31 = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41
    b32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42
    b33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41
    b34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32
    b41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42
    b42 = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41
    b43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42
    b44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31
    imatrix = 1./detA*np.array([[b11,b12,b13,b14],[b21,b22,b23,b24],[b31,b32,b33,b34],[b41,b42,b43,b44]])
    return imatrix,detA

def calc_M(a,rho):
    ai=a[:,0]; aq=a[:,1];au=a[:,2];av=a[:,3] #= a# = a[1]; au = a[2]; av = a[3]
#    onopol=np.exp(-ai*x)
    rhoq = rho[:,0]; rhou = rho[:,1]; rhov = rho[:,2]
    a2 = aq**2.+au**2.+av**2.
    p2 = rhoq**2.+rhou**2.+rhov**2.
    if np.sum(a2)==0. and np.sum(p2)==0.:
        return np.identity(4)*onopol,0.,0.,0.,0.
    else:
        ap = aq*rhoq+au*rhou+av*rhov
        lam1 = np.sqrt(np.sqrt((a2-p2)**2./4.+ap**2.)+(a2-p2)/2.)
        lam2 = np.sqrt(np.sqrt((a2-p2)**2./4.+ap**2.)-(a2-p2)/2.)
        theta = lam1**2.+lam2**2.
        sig = np.sign(ap)
        zero = 0.*ap
        M2 = np.array([[zero,lam2*aq-sig*lam1*rhoq,lam2*au-sig*lam1*rhou,lam2*av-sig*lam1*rhov],[lam2*aq-sig*lam1*rhoq,zero,sig*lam1*av+lam2*rhov,-sig*lam1*au-lam2*rhou],[lam2*au-sig*lam1*rhoq,-sig*lam1*av-lam2*rhov,zero,sig*lam1*aq+lam2*rhoq],[lam2*av-sig*lam1*rhov,sig*lam1*au+lam2*rhou,-sig*lam1*aq-lam2*rhoq,zero]])
        M3 = np.array([[zero,lam1*aq+sig*lam2*rhoq,lam1*au+sig*lam2*rhou,lam1*av+sig*lam2*rhov],[lam1*aq+sig*lam2*rhoq,zero,-sig*lam2*av+lam1*rhov,sig*lam2*au-lam1*rhou],[lam1*au+sig*lam2*rhou,sig*lam2*av-lam1*rhov,zero,-sig*lam2*aq+lam1*rhoq],[lam1*av+sig*lam2*rhov,-sig*lam2*au+lam1*rhou,sig*lam2*aq-lam1*rhoq,zero]])
        M4 = np.array([[(a2+p2)/2.,av*rhou-au*rhov,aq*rhov-av*rhoq,au*rhoq-aq*rhou],[au*rhov-av*rhou,aq**2.+rhoq**2.-(a2+p2)/2.,aq*au+rhoq*rhou,av*aq+rhov*rhoq],[av*rhoq-aq*rhov,aq*au+rhoq*rhou,au**2.+rhou**2.-(a2+p2)/2.,au*av+rhou*rhov],[aq*rhou-au*rhoq,av*aq+rhoq*rhov,av*au+rhou*rhov,av**2.+rhov**2.-(a2+p2)/2.]])
    return M2,M3,M4,lam1,lam2,theta

# this is intended to be called with 4x4 matrices M1,M2,M3,M4 and scalar lam1,lam2,x
def calc_O_M(M1,M2,M3,M4,lam1,lam2,ai,x):
    onopol = np.exp(-ai*x)
    O = onopol*(1./2.*(np.cosh(lam1*x)+np.cos(lam2*x))*M1 - np.sin(lam2*x)*M2-np.sinh(lam1*x)*M3+1./2.*(np.cosh(lam1*x)-np.cos(lam2*x))*M4)
    return O

def calc_O(a,rho,x):
#    onopol = np.exp(-a[0]*x)
    ai,aq,au,av = a# = a[1]; au = a[2]; av = a[3]
    onopol = np.exp(-ai*x)
    rhoq,rhou,rhov = rho#rho[0]; rhou = rho[1]; rhov = rho[2]
    a2 = aq**2.+au**2.+av**2.
    p2 = rhoq**2.+rhou**2.+rhov**2.
    if np.sum(a2)==0. and np.sum(p2)==0.:
        return np.identity(4)*onopol,0.,0.,0.,0.
    else:
        ap = aq*rhoq+au*rhou+av*rhov
        lam1 = np.sqrt(np.sqrt((a2-p2)**2./4.+ap**2.)+(a2-p2)/2.)
        lam2 = np.sqrt(np.sqrt((a2-p2)**2./4.+ap**2.)-(a2-p2)/2.)
        theta = lam1**2.+lam2**2.
        sig = np.sign(ap)
        M1 = np.identity(4)
        M2 = 1./theta*np.array([[0.,lam2*aq-sig*lam1*rhoq,lam2*au-sig*lam1*rhou,lam2*av-sig*lam1*rhov],[lam2*aq-sig*lam1*rhoq,0.,sig*lam1*av+lam2*rhov,-sig*lam1*au-lam2*rhou],[lam2*au-sig*lam1*rhoq,-sig*lam1*av-lam2*rhov,0.,sig*lam1*aq+lam2*rhoq],[lam2*av-sig*lam1*rhov,sig*lam1*au+lam2*rhou,-sig*lam1*aq-lam2*rhoq,0.]])
        M3 = 1./theta*np.array([[0.,lam1*aq+sig*lam2*rhoq,lam1*au+sig*lam2*rhou,lam1*av+sig*lam2*rhov],[lam1*aq+sig*lam2*rhoq,0.,-sig*lam2*av+lam1*rhov,sig*lam2*au-lam1*rhou],[lam1*au+sig*lam2*rhou,sig*lam2*av-lam1*rhov,0.,-sig*lam2*aq+lam1*rhoq],[lam1*av+sig*lam2*rhov,-sig*lam2*au+lam1*rhou,sig*lam2*aq-lam1*rhoq,0.]])
        M4 = 2./theta*np.array([[(a2+p2)/2.,av*rhou-au*rhov,aq*rhov-av*rhoq,au*rhoq-aq*rhou],[au*rhov-av*rhou,aq**2.+rhoq**2.-(a2+p2)/2.,aq*au+rhoq*rhou,av*aq+rhov*rhoq],[av*rhoq-aq*rhov,aq*au+rhoq*rhou,au**2.+rhou**2.-(a2+p2)/2.,au*av+rhou*rhov],[aq*rhou-au*rhoq,av*aq+rhoq*rhov,av*au+rhou*rhov,av**2.+rhov**2.-(a2+p2)/2.]])
        O = onopol*(1./2.*(np.cosh(lam1*x)+np.cos(lam2*x))*M1 - np.sin(lam2*x)*M2-np.sinh(lam1*x)*M3+1./2.*(np.cosh(lam1*x)-np.cos(lam2*x))*M4)
        return O,M1,M2,M3,M4

def num(x=np.array([1.]),j=np.array([1.,0.,0.,0.]),a=np.array([1.,0.,0.,0.]),rho=np.array([0.,0.,0.]),I0=np.array([0.,0.,0.,0.])):
    integrand = np.zeros((len(x),4))
    if len(x) < 2:
        dx = x
    else:
        dx = np.append(x[0],x[1:]-x[0:-1])
    if len(np.shape(a)) < 2:
# reform a,rho,j arrays to be of right size
        a = np.tile(a,len(x)).reshape(len(x),4)
        rho = np.tile(rho,len(x)).reshape(len(x),3)
        j = np.tile(j,len(x)).reshape(len(x),4)
        
    i = np.zeros((len(x),4))
    intprev = np.zeros(4); iprev = I0; xprev = 0.; jprev = np.zeros(4)
    i[0,:] = I0
    for k in range(len(x)-1):
        K = opacity_matrix(a[k,:],rho[k,:])
        K1 = opacity_matrix(a[k+1,:],rho[k+1,:])
        dIds = j[k,:]-K.dot(iprev)
# "symplectic" attempt:
        inew=dIds*(x[k+1]-x[k])+i[k,:]
        dIds1 = j[k,:]-K.dot(inew)
        dIds[1]=dIds1[1]
        i[k+1,:] = dIds*(x[k+1]-x[k])+i[k,:]
        iprev = i[k+1,:]
        integrand[k,:] = dIds

    return i,integrand


# calculate intensity over some set of coefficients j,a,rho at positions x for initial intensity I0
def intensity(x=np.array([1.]),j=np.array([1.,0.,0.,0.]),a=np.array([1.,0.,0.,0.]),rho=np.array([0.,0.,0.]),I0=np.array([0.,0.,0.,0.])):
    o = np.zeros((len(x),4,4))
    integrand = np.zeros((len(x),4))
    if len(x) < 2:
        dx = x
    else:
        dx = np.append(x[0],x[1:]-x[0:-1])
    if len(np.shape(a)) < 2:
# reform a,rho,j arrays to be of right size
        a = np.tile(a,len(x)).reshape(len(x),4)
        rho = np.tile(rho,len(x)).reshape(len(x),3)
        j = np.tile(j,len(x)).reshape(len(x),4)
        
    i = np.zeros((len(x),4))
    intprev = np.zeros(4); iprev = I0; xprev = 0.; jprev = np.zeros(4)
#    xx = np.append(0.,x)
#    for k in range(len(x)-1):
#        oxk,M1,M2,M3,M4 = calc_O(a[k,:],rho[k,:],x[-1]-x[k])
#        oxk1,M1,M2,M3,M4 = calc_O(a[k+1],rho[k+1,:],x[-1]-x[k+1])
#        o[k,:,:],M1,M2,M3,M4 = calc_O(a[k,:],rho[k,:],xx[-k-1])
#        integrand[k,:] = o[k,:,:].dot(jprev)
#        i1 = oxk.dot(j[k,:])
#        i2 = oxk1.dot(j[k+1,:])
#        i[k+1,:] = (i1+i2)/2.*(x[k+1]-x[k])+oxk.dot(i[k,:])
#        iprev = i[k,:]; xprev = x[k]; intprev = integrand[k,:]; jprev = j[k,:]
        
# intensity for constant coefs is integral along path + attenuated initial intensity
    for k in range(len(x)):
        o[0,:,:],M1,M2,M3,M4 = calc_O(a[k,:],rho[k,:],1.)
        integrand[k,:] = o[k,:,:].dot(j[k,:])

#    i = np.append(np.zeros((1,4)),scipy.integrate.cumtrapz(integrand,np.transpose(np.tile(x,4).reshape((4,len(x)))),axis=0),axis=0) + o[0,:,:].dot(I0)
    intatten = o[0,:,:].dot(I0)
    for k in range(4):
        i[:,k] = np.append(0.,scipy.integrate.cumtrapz(integrand[:,k],x)) + intatten[k]

    return i,o,dx,integrand

# intensity over some set of coefficients j,a,rho at positions x for initial intensity I0 for arbitrary coefficients
def intensity_var(x=np.array([1.]),j=np.array([1.,0.,0.,0.]),a=np.array([1.,0.,0.,0.]),rho=np.array([0.,0.,0.]),I0=np.array([0.,0.,0.,0.])):
    o = np.zeros((len(x),4,4))
    integrand = np.zeros((len(x),4))
    if len(x) < 2:
        dx = x
    else:
        dx = np.append(x[0],x[1:]-x[0:-1])
    if len(np.shape(a)) < 2:
# reform a,rho,j arrays to be of right size
        a = np.tile(a,len(x)).reshape(len(x),4)
        rho = np.tile(rho,len(x)).reshape(len(x),3)
        j = np.tile(j,len(x)).reshape(len(x),4)
        
    i = np.zeros((len(x),4))
    intprev = I0; iprev = I0; xprev = 0.; jprev = np.zeros(4)
# intensity for constant coefs is integral along path + attenuated initial intensity
    identity = np.identity(4)
    o[0,:,:] = identity
    M1=identity
    print('calc_M: ',np.shape(a),np.shape(rho))
    M2,M3,M4,lam1,lam2,theta = calc_M(a,rho)
    for k in range(len(x)-1):
#        o[k+1,:,:],M1,M2,M3,M4 = calc_O(a[k,:],rho[k,:],x[k+1]-x[k])
        o[k+1,:,:] = calc_O_M(M1,1./theta[k]*M2[:,:,k],1./theta[k]*M3[:,:,k],2./theta[k]*M4[:,:,k],lam1[k],lam2[k],a[k,0],x[k+1]-x[k])
        jj=j[k,:]
        oo=o[k,:,:]
# try "symplectic" where intprev is updated for Q early:
#        iupdate = oo.dot(jj*(x[k+1]-x[k])+intprev)
#        intprev[1]=iupdate[1]
        i[k+1,:] = oo.dot(jj)*(x[k+1]-x[k])+oo.dot(intprev)
        intprev = i[k+1,:]

    return i,o,dx

def intensity_var_backwards(x=np.array([1.]),j=np.array([1.,0.,0.,0.]),a=np.array([1.,0.,0.,0.]),rho=np.array([0.,0.,0.]),I0=np.array([0.,0.,0.,0.])):
    o = np.zeros((len(x),4,4))
    ocum = np.zeros((len(x),4,4))
    integrand = np.zeros((len(x),4))
    if len(x) < 2:
        dx = x
    else:
        dx = np.append(x[0],x[1:]-x[0:-1])
    if len(np.shape(a)) < 2:
# reform a,rho,j arrays to be of right size
        a = np.tile(a,len(x)).reshape(len(x),4)
        rho = np.tile(rho,len(x)).reshape(len(x),3)
        j = np.tile(j,len(x)).reshape(len(x),4)
        
    i = np.zeros((len(x),4))
    intcur = np.zeros(4)
# intensity for constant coefs is integral along path + attenuated initial intensity
    identity = np.identity(4)
#    o[0,:,:],M1,M2,M3,M4 = calc_O(a[0,:],rho[0,:],x[1]-x[0])
    ocum[0,:,:]=np.identity(4)
    for k in np.arange(len(x)-1)+1:
        o[k,:,:],M1,M2,M3,M4 = calc_O(a[k-1,:],rho[k-1,:],x[k]-x[k-1])
# JAD 3/22 is this really a dot? I guess it is a matrix multiplication
#        ocum[k,:,:]=ocum[k-1,:,:].dot(o[k,:,:])
        ocum[k,:,:]=np.matmul(ocum[k-1,:,:],o[k,:,:],out=ocum[k,:,:])
        jj = j[k,:]
        integrand[k,:] = ocum[k,:,:].dot(jj)

    print('len: ',len(i[:,0]),len(integrand[:,0]), len(x))

    for m in range(4):
#        i[:,m] = scipy.integrate.cumtrapz(integrand[:,m],x,initial=0.)
        i[:,m] = np.cumsum(integrand[:,m]*dx)

    return i,o,dx,integrand


def invert_delo_matrix_thin(dx,K,ki,delta):
#    matrix = np.identity(4)*(1.-delta/2.)+0.5*dx*K
    matrix = np.identity(4)*(1.-delta/2.+delta**2./6.)+(0.5*dx-1./6.*dx**2.*ki)*K
    imatrix = np.linalg.inv(matrix)
    return matrix,imatrix

def calc_delo_P_thin(imatrix,dx,j,j1,ki,ki1):
#    return imatrix.dot(0.5*dx*j+0.5*j1*dx)
    return imatrix.dot((0.5*dx*j-1./6.*dx**2.*ki*j)+(0.5*j1*dx-1./3.*dx**2.*ki*j1))
#    return imatrix.dot(dx*j)

def calc_delo_Q_thin(imatrix,dx,ki,ki1,K1):
#    return imatrix.dot(np.identity(4)*(1.-0.5*dx*ki)-0.5*ki/ki1*dx*K1)
#    return imatrix.dot(np.identity(4)*(1.-0.5*dx*ki)-0.5*dx*K1)
    return imatrix.dot(np.identity(4)*(1.-0.5*dx*ki+1./6.*dx**2.*ki**2.)-(0.5*dx-1./3.*dx**2.)*K1)
#    return np.identity(4)

def invert_delo_matrix(F,G,Kp):
    matrix = np.identity(4)+(F-G)*Kp
    imatrix = np.linalg.inv(matrix)
    return matrix,imatrix

def calc_delo_P(imatrix,F,G,Sp,Sp1):
    return imatrix.dot(((F-G)*Sp+G*Sp1))

def calc_delo_Q(imatrix,E,F,G,Kp1):
    return imatrix.dot(np.identity(4)*E-G*Kp1)

def delo_intensity(dx=np.array([1.]),j=np.array([1.,0.,0.,0.]),a=np.array([1.,0.,0.,0.]),rho=np.array([0.,0.,0.]),I0=np.array([0.,0.,0.,0.]),thin=1e-2):
    x = np.append(0.,np.cumsum(dx))
#    if len(dx) < 2:
#        dx = x
#    else:
#        dx = np.append(x[0],x[1:]-x[0:-1])
    if len(np.shape(a)) < 2:
# reform a,rho,j arrays to be of right size
        a = np.tile(a,len(x)).reshape(len(x),4)
        rho = np.tile(rho,len(x)).reshape(len(x),3)
        j = np.tile(j,len(x)).reshape(len(x),4)
        
#    ki4 = np.transpose(np.tile(a[:,0],4).reshape(4,len(x)))
    i = np.zeros((len(x),4)); Q = np.zeros((len(x),4,4)); P = np.zeros((len(x),4)); im = np.zeros((len(x),4,4))
    QQ = np.zeros((len(x),4,4)); PP = np.zeros((len(x),4)); imm = np.zeros((len(x),4,4)); ii = np.zeros((len(x),4))
#    i[0,:] = I0
#    tau = -(scipy.integrate.cumtrapz(x[-1]-x,a[::-1,0]))[::-1]
    tau = np.append(0.,-scipy.integrate.cumtrapz(a[:,0],x[::-1]))
#    delta = a[:,0]*dx
    delta = tau[1:] - tau[0:-1]
    E = np.exp(-delta)
    F = 1.-E
    G = (1.-(1.+delta)*E)/delta
# opt thin version to avoid errors from 1/delta w/ delta --> 0
    Gt = 0.5*delta
    Ft = delta
    Et = 1.-delta
# integration is from deepest point out for starting intensity I0
    i[-1,:] = I0; ii[-1,:] = I0; iprev = I0; iprevt = I0
    for k in (range(len(x)-1))[::-1]:
#        print 'k: ',k,len(F),len(G),len(delta)
        K = opacity_matrix(a[k,:],rho[k,:])
        K1 = opacity_matrix(a[k+1,:],rho[k+1,:])
        Sp = j[k,:]/a[k,0]
        Sp1 = j[k+1,:]/a[k+1,0]
        Kp = K/a[k,0]-np.identity(4); Kp1 = K1/a[k+1,0]-np.identity(4)
        matrix,imatrix = invert_delo_matrix(F[k],G[k],Kp)
        mt,imt = invert_delo_matrix_thin(dx[k],K,a[k,0],delta[k])
        pt = calc_delo_P(imatrix,F[k],G[k],Sp,Sp1)
        qt = calc_delo_Q(imatrix,E[k],F[k],G[k],Kp1)
        ptt = calc_delo_P_thin(imt,dx[k],j[k,:],j[k+1,:],a[k,0],a[k+1,0])
        qtt = calc_delo_Q_thin(imt,dx[k],a[k,0],a[k+1,0],K1)
        mtt,imtt = invert_delo_matrix(Ft[k],Gt[k],Kp)
#        pttt = calc_delo_P(imatrix,Ft[k],Gt[k],Sp,Sp1)
#        qttt = calc_delo_Q(imatrix,Et[k],Ft[k],Gt[k],Kp1)
        if delta[k] > thin:
            i[k,:] = pt + qt.dot(iprev)
        else:
            i[k,:] = ptt + qtt.dot(iprev)
        P[k,:] = pt
        Q[k,:,:] = qt
        im[k,:,:] = imatrix
        imm[k,:,:] = imt
        QQ[k,:,:] = qtt
        PP[k,:] = ptt
        iprev = i[k,:]
        iprevt = ii[k,:]
        
    return i,Q,P,im,delta,dx,QQ,PP,imm,ii
