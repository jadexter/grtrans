import constants as pc
import numpy as np

N1 = 288
N2 = 128
N3 = 1
#N1 = 18
#N2 = 8
#N3 = 8
deltaphi = 2.*pc.pi
R0 = 0.3
Rout = 1e5

fracdisk = 0.36
fracjet = 0.3
fraccorona = 1. - fracdisk - fracjet
nu = 0.75
rsjet = 0.
a = 0.99
rhor = 1.+np.sqrt(1.-a*a)
Rin = 0.83*rhor
r0grid = 5.*Rin
r0jet = Rin
rjetend = 5.
r0disk = Rin
hoverr = 0.1
rmaxdisk=34.
rdiskend = 300.
logfac = 0.
rinsidediskmax = (1.+logfac*np.log10(rdiskend/rmaxdisk))**(2./nu)
rbr = 1000.
npow2 = 4.
cpow2 = 1.
x10 = 2.4
x20 = -1.+1./N2
x1min = np.log(Rin-R0)
x1br = np.log(rbr-R0)

def step(x):
    return 0.5*(np.sign(x)+1.)

def fcore(x):
    return 0.5*((x+1.)+2./128./pc.pi*(-140.*np.sin(pc.pi*(x+1.)/2.)+10./3.*np.sin(3.*pc.pi*(x+1.)/2.)+2./5.*np.sin(5.*pc.pi*(x+1.)/2.)))

def f(x):
#    x = np.array([xx])
#    zero = np.zeros(np.size(x))
    result = np.where(x > -1.,np.where(x < 1.,fcore(x),x),0.)
    return result

def limlin(x,x0,dx):
    return x0-dx*f(-(x-x0)/dx)

def minlin(x,x0,dx,y0):
    return y0+dx*f((x-x0)/dx)

def mins(f1,f2,df):
    return limlin(f1,f2,df)

def maxs(f1,f2,df):
    return -limlin(-f1,-f2,df)

def minmaxs(f1,f2,df,dir):
    s = step(dir)
    return s*maxs(f1,f2,df)+(1.-s)*mins(f1,f2,df)

def calc_r(x1):
    return R0+np.exp(x1+cpow2*(step(x1-x1br)*(x1-x1br))**npow2)

def rinsidedisk(x1):
    r = calc_r(x1)
    return (1.+logfac*np.log10(mins(maxs(1.,r/rmaxdisk,0.5),rdiskend/rmaxdisk,0.5*rdiskend/rmaxdisk)))**(2./nu)

def rbeforedisk(x1):
    return mins(calc_r(x1),r0disk,0.5*r0disk)

def rafterdisk(x1):
    return maxs(1.,1.+(calc_r(x1)-rdiskend)*r0jet/(rjetend*r0disk*rinsidediskmax),0.5*rdiskend*r0jet/(rjetend*r0disk*rinsidediskmax))

def fakerdisk(x1,x2):
    return rbeforedisk(x1)*rinsidedisk(x1)*rafterdisk(x1)

def fakerjet(x1,x2):
    r = calc_r(x1)
    return mins(r,r0jet,0.5*r0jet)*maxs(1.,r/rjetend,0.5)

def Ftr1(x):
    return 1./128.*(64.+np.cos(5./2.*pc.pi*(1.+x))+70.*np.sin(pc.pi*x/2.)+5.*np.sin(3.*pc.pi*x/2.))

def Ftr(x):
#    x = np.array([xx])
#    zero = np.zeros(np.size(x))
#    one = np.zeros(np.size(x))
    return np.where(x < 0.,0.,np.where(x > 1.,1.,Ftr1(2.*x-1.)))

def f0(x,xa,xb,r0a,r0b):
    return r0a+(r0b-r0a)*Ftr((x-xa)/(xb-xa))

def fac(x2):
    return f0(np.abs(x2),fracdisk,1.-fracjet,0.,1.)

def faker(x1,x2):
    return fakerdisk(x1,x2)*(1.-fac(x2))+fakerjet(x1,x2)*fac(x2)-rsjet*Rin

def prefact(x1,x2):
    return (faker(x1,x2)/r0grid)**(nu/2.)

#def thori(x1,x2):
#    th = pc.pi/2.+np.arctan(np.tan(x2*pc.pi/2.)*prefact(x1,x2))
#    return th

def prefactdisk(x1,x2):
    return (fakerdisk(x1,x2)/r0grid)**(nu/2.)

def prefactjet(x1,x2):
    return (fakerjet(x1,x2)/r0grid)**(nu/2.)

def thdisk(x1,x2):
    return pc.pi/2.+np.arctan(np.tan(x2*pc.pi/2.)*prefactdisk(x1,x2))

def thjet(x1,x2):
    return pc.pi/2.+np.arctan(np.tan(x2*pc.pi/2.)*prefactjet(x1,x2))

def thori(x1,x2):
    return thdisk(x1,x2)*(1.-fac(x2))+thjet(x1,x2)*fac(x2)-rsjet*Rin

def th0(x1,x2):
    return np.arcsin(calc_r(x10)*np.sin(thori(x10,x20))/calc_r(x1))

def th1in(x1,x2):
    return np.arcsin(calc_r(x10)*np.sin(thori(x10,x2))/calc_r(x1))

def th2in(x1,x2):
    th012 = th0(x1,x2)
    return ((thori(x1,x2)-thori(x1,x20))/(thori(x1,0.)-thori(x1,x20)))*(thori(x1,0.)-th012)+th012

def func1(x1,x2):
    return np.sin(thori(x1,x2))

def func2(x1,x2):
    arg1 = np.real(np.sin(th1in(x1,x2)))
    arg2 = np.real(np.sin(th2in(x1,x2)))
    arg3 = np.abs(np.real(np.sin(th2in(x1,-1.)))-np.real(np.sin(th1in(x1,-1.))))
    arg4 = x1 - x10
    return minmaxs(arg1,arg2,arg3,arg4)

def dfunc(x1,x2):
    return func2(x1,x2)-func1(x1,x2)

def calc_th(x1,x2):
    r = calc_r(x1)
    argmax = calc_r(np.log(0.5*np.exp(x10)+0.5*np.exp(x1min)))*(np.abs(dfunc(np.log(0.5*np.exp(x10)+0.5*np.exp(x1min)),x2)))
    th = np.arcsin(maxs(r*func1(x1,x2),r*func2(x1,x2),argmax)/r)
    return th

def x1func(i,N,x1min,x1max):
    return x1min+(x1max-x1min)*(i+0.5)/N

def construct_grid():
# r grid
    x1 = np.zeros((N1,N2,N3),dtype=complex)
    x2 = np.zeros((N1,N2,N3),dtype=complex)
    x3 = np.zeros((N1,N2,N3),dtype=complex)
    i = np.arange(N1)
    x1max = 8.25133
#    x1max = findx1(Rout)
    x1d = x1func(i,N1,x1min,x1max)
    x2d = (np.arange(N2)+0.5)/N2
    x3d = (np.arange(N3)+0.5)/N3
    for i in range(N1):
        for j in range(N2):
            x3[i,j,:] = x3d
    for i in range(N1):
        for j in range(N3):
            x2[i,:,j] = x2d
    for i in range(N2):
        for j in range(N3):
            x1[:,i,j]=x1d
    
    r = calc_r(x1.flatten())
    th = calc_th(x1.flatten(),x2.flatten())
    phi = x3*2.*pc.pi
    return x1,x2,x3,r,th,phi
