from radtrans_integrate import radtrans_integrate
import scipy.integrate
import numpy as np

# analytic solution to pol rad trans equations with emission / absorption only in Stokes I,Q
def analytic_radtrans_absorption(ji,jq,ai,aq,s):
    alpha = ai+aq
    denom=1./2./(ai-aq)/(ai+aq)
    num=(ji*aq-jq*ai)*(np.exp(-alpha*s)-np.exp(2.*aq*s-alpha*s))+(jq*aq-ji*ai)*(np.exp(-alpha*s)+np.exp(2.*aq*s-alpha*s)-2.)
    ii=num*denom
    numq=((jq*aq-ji*ai)*(np.exp(-alpha*s)-np.exp(2.*aq*s-alpha*s))+(ji*aq-jq*ai)*(np.exp(-alpha*s)+np.exp(-alpha*s+2.*aq*s)-2.))
    iq=numq*denom
    return ii,iq

# analytic solution to pol rad trans equations with emission in Q, Faraday rotation and conversion
def analytic_radtrans_faraday(jq,ju,jv,rhoq,rhov,s):
# redone from notes 4/13/2015 solving ODEs with all pol emissivities, just assuming rhou = 0
    rho=np.sqrt(rhoq**2.+rhov**2.)
#    iq=rhov**2.*(rhoq+rhov)/rho**4.*ju*(np.cos(rho*s)-1.)+rhoq*rhov/rho**3.*(jq-jv)*np.sin(rho*s)+rhoq*s/rho**2.*(jq*rhoq+jv*rhov)
    iq = ju*rhov/rho**2.*(np.cos(rho*s)-1.)-rhov/rho**3.*(jv*rhoq-jq*rhov)*np.sin(rho*s)+rhoq*s/rho**2.*(jq*rhoq+jv*rhov)
#    iu=-rhov*(rhoq+rhov)/rho**3.*ju*np.sin(rho*s)+rhoq/rho**2.*(jq-jv)*(1.-np.cos(rho*s))
    iu=ju/rho*np.sin(rho*s)+(jq*rhov-jv*rhoq)/rho**2.*(1.-np.cos(rho*s))
#    iv=rhov*s/rho**2.*(jq*rhoq+jv*rhov)+rhoq*rhov*(rhoq+rhov)/rho**4.*ju*(1.-np.cos(rho*s))+rhoq**2./rho**3.*(jq-jv)*np.sin(rho*s)
    iv=(1.-np.cos(rho*s))*ju*rhoq/rho**2.-np.sin(rho*s)*(jq*rhov-jv*rhoq)*rhoq/rho**3.+rhov*s*(jq*rhoq+jv*rhov)/rho**2.
    return iq,iu,iv

# full solution to pol rad trans equations with emission in QUV, Faraday rotation and conversion
# see notes or mathematica analytic_radtrans_polabs_fcfr.nb
def analytic_radtrans_faraday_full(jq,ju,jv,rhoq,rhout,rhov,s):
# account for sign error made in equation 47 of Dexter 2016 for rhoU
    rhou=-rhout
    rho=np.sqrt(rhoq**2.+rhou**2.+rhov**2.)
    crho=np.cos(rho*s); srho=np.sin(rho*s)
    iq=(crho-1.)*(jv*rhou+ju*rhov)/rho**2.+srho/rho**3.*(jq*(rhou**2.+rhov**2.)-jv*rhoq*rhov+ju*rhoq*rhou)+rhoq*s/rho**2.*(jq*rhoq-ju*rhou+jv*rhov)
    iu=(crho-1.)*(jv*rhoq-jq*rhov)/rho**2.+srho/rho**3.*(ju*(rhoq**2.+rhov**2.)+jv*rhou*rhov+jq*rhoq*rhou)+rhou*s/rho**2.*(ju*rhou-jq*rhoq-jv*rhov)
    iv=(1.-crho)*(ju*rhoq+jq*rhou)/rho**2.+srho/rho**3.*(jv*(rhou**2.+rhoq**2.)+ju*rhou*rhov-jq*rhoq*rhov)+rhov*s/rho**2.*(jv*rhov+jq*rhoq-ju*rhou)
    return iq,iu,iv

def run(abs=1,far=1,bfar=1,save=0):
    # do strong FR test with emission and compare to analytic solution for LSODA, DELO, FORMAL methods
    x = np.cumsum(np.zeros(3000)+1e-2)-1e-2
    rhoarr = np.zeros((3000,3))
    rhoarr[:,2]=10.#; rhoarr[:,1]=-4.
    jarr = np.tile(np.array([1.,0.7,0.,0.]),3000).reshape(3000,4)
#a = np.array([10.,4.,1.,7.])
#a2 = np.array([5.,3.,2.,1.])
#a2arr = np.tilbe(a2,150).reshape(150,4)
    aarr = np.zeros((3000,4))
#aarr = np.append(a2arr/3.,3.*aarr,axis=0)
    Karr = (np.append(aarr,rhoarr,axis=1))
    tau = np.append(0.,scipy.integrate.cumtrapz(Karr[:,0],x))
    radtrans_integrate.init_radtrans_integrate_data(0,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ilb = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(1,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    idb = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(2,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ifb = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
# analytic solution
    ii = jarr[:,0]*x
    iu = jarr[:,1]/rhoarr[:,2]*(1.-np.cos(x*rhoarr[:,2]))
    iq = jarr[:,1]/rhoarr[:,2]*np.sin(x*rhoarr[:,2])
    print np.shape(iu), np.shape(jarr), np.shape(ifb)
    print 'max error faraday basic lsoda: ',np.max(np.abs(ii-ilb[0,:])),np.max(np.abs(iq-ilb[1,:])),np.max(np.abs(iu-ilb[2,:]))
    print 'max error faraday basic formal: ',np.max(np.abs(ii-ifb[0,:])),np.max(np.abs(iq-ifb[1,:])),np.max(np.abs(iu-ifb[2,:]))
    print 'max error faraday basic delo: ',np.max(np.abs(ii-idb[0,:])),np.max(np.abs(iq-idb[1,:])),np.max(np.abs(iu-idb[2,:]))

# now with both Faraday rotation and conversion
    rhoarr[:,0] = -4.; rhoarr[:,2]=10.
    aarr[:,:] = 0.; jarr[:,:]=0.
    jarr[:,0]=1e-8; jarr[:,1]=0.1; jarr[:,3]=0.1; jarr[:,2]=0.1
    Karr = (np.append(aarr,rhoarr,axis=1))
    radtrans_integrate.init_radtrans_integrate_data(0,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ilf = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(1,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    idf = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(2,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    iff = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
# analytic solution
    ii = jarr[:,0]*x
    iq,iu,iv = analytic_radtrans_faraday(jarr[:,1],jarr[:,2],jarr[:,3],rhoarr[:,0],rhoarr[:,2],x)
    print np.shape(iu), np.shape(jarr), np.shape(iff)
    print 'max error faraday lsoda: ',np.max(np.abs(ii-ilf[0,:])),np.max(np.abs(iq-ilf[1,:])),np.max(np.abs(iu-ilf[2,:])),np.max(np.abs(iv-ilf[3,:]))
    print 'max error faraday formal: ',np.max(np.abs(ii-iff[0,:])),np.max(np.abs(iq-iff[1,:])),np.max(np.abs(iu-iff[2,:])),np.max(np.abs(iv-iff[3,:]))
    print 'max error faraday delo: ',np.max(np.abs(ii-idf[0,:])),np.max(np.abs(iq-idf[1,:])),np.max(np.abs(iu-idf[2,:])),np.max(np.abs(iv-idf[3,:]))

# with rhou as well
    rhoarr[:,0] = -4.; rhoarr[:,2]=100.; rhoarr[:,1]=7.
    rhoarr*=1000.
    aarr[:,:] = 0.; jarr[:,:]=0.
    jarr[:,0]=1e-8; jarr[:,1]=0.1; jarr[:,3]=0.1; jarr[:,2]=0.1
    Karr = (np.append(aarr,rhoarr,axis=1))
    radtrans_integrate.init_radtrans_integrate_data(0,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ilf2 = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(1,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    idf2 = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(2,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    iff2 = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
# analytic solution
    ii = jarr[:,0]*x
    iq,iu,iv = analytic_radtrans_faraday_full(jarr[:,1],jarr[:,2],jarr[:,3],rhoarr[:,0],rhoarr[:,1],rhoarr[:,2],x)
    print np.shape(iu), np.shape(jarr), np.shape(iff2)
    print 'max error faraday full lsoda: ',np.max(np.abs(ii-ilf2[0,:])),np.max(np.abs(iq-ilf2[1,:])),np.max(np.abs(iu-ilf2[2,:])),np.max(np.abs(iv-ilf2[3,:]))
    print 'max error faraday full formal: ',np.max(np.abs(ii-iff2[0,:])),np.max(np.abs(iq-iff2[1,:])),np.max(np.abs(iu-iff2[2,:])),np.max(np.abs(iv-iff2[3,:]))
    print 'max error faraday full delo: ',np.max(np.abs(ii-idf2[0,:])),np.max(np.abs(iq-idf2[1,:])),np.max(np.abs(iu-idf2[2,:])),np.max(np.abs(iv-idf2[3,:]))

    
# now do strong absorption test with emission the same way
    aarr[:,0] = 5.; aarr[:,1] = 4.; jarr[:,:] = 0.
    jarr[:,0] = 1.; jarr[:,1] = 0.8
    rhoarr[:,:] = 0.
    Karr = (np.append(aarr,rhoarr,axis=1))
    tau = np.append(0.,scipy.integrate.cumtrapz(Karr[:,0],x))
    radtrans_integrate.init_radtrans_integrate_data(0,4,3000,3000,2000.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ila = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(1,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ida = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
    radtrans_integrate.init_radtrans_integrate_data(2,4,3000,3000,10.,0.1,1e-8,1e-6,1e-2)
    radtrans_integrate.integrate(x[::-1],jarr[:,:],Karr[:,:],tau,4)
    ifa = radtrans_integrate.intensity.copy()
    radtrans_integrate.del_radtrans_integrate_data()
# analytic solution
    ii,iqq = analytic_radtrans_absorption(jarr[:,0],jarr[:,1],aarr[:,0],aarr[:,1],x)
#iu = jarr[:,1]/rhoarr[:,2]*(1.-np.cos(x*rhoarr[:,2]))
#iq = jarr[:,1]/rhoarr[:,2]*np.sin(x*rhoarr[:,2])
    print np.shape(iu), np.shape(jarr), np.shape(ifa)
    print 'max error absorption lsoda: ',np.max(np.abs(ii-ila[0,:])),np.max(np.abs(iqq-ila[1,:]))
    print 'max error absorption formal: ',np.max(np.abs(ii-ifa[0,:])),np.max(np.abs(iqq-ifa[1,:]))
    print 'max error absorption delo: ',np.max(np.abs(ii-ida[0,:])),np.max(np.abs(iqq-ida[1,:]))

    if save==1:
        analytic = np.array([iq,iu,iv,ii,iqq,x])
        num = np.append(np.append(np.append(np.append(np.append(ilf,idf),iff),ila),ida),ifa)
        np.savetxt('unit_tests_integration_output_analytic.txt',analytic)
        np.savetxt('unit_tests_integration_output_num.txt',num)

    return ilb,idb,ifb,ilf,idf,iff,ila,ida,ifa,ilf2,idf2,iff2,iq,iu,iv
