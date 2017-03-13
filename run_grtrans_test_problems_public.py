import grtrans_batch as gr
import pickle
import numpy as np
import copy

def run_test_problems(save=0,pgrtrans=0,nosphacc=0,compile=0):
    # run grtrans test problems
    tol=1e-2; failed=[]; xlist=[]
    xlist.append(gr.grtrans())
    if compile > 0:
        if pgrtrans==0:
            xlist[-1].compile_grtrans()
        else:
            xlist[-1].compile_pgrtrans()
    passed=0; max_passed=0
    if nosphacc <= 0:
# sphacc
#    xlist[-1].write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=15,nmu=1,fmin=2.41e10,fmax=6.31e14,ename='POLSYNCHTH',nvals=4,spin=0.,mbh=10.,standard=1,uout=.003,nn=[100,100,100")
# New tests of 1d intensity profile & full spectrum 12/14/2012
        xlist[-1].write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=25,nmu=1,fmin=1e8,fmax=1e15,ename='SYNCHTHAV',nvals=1,spin=0.,mbh=1.,standard=1,nn=[10000,1,100],gridvals=[0.,400.,0.,0.],uout=.0025,oname='sphacc_abs.out')
        if pgrtrans==0:
            xlist[-1].run_grtrans()
            xlist[-1].read_grtrans_output()
        else:
            xlist[-1].run_pgrtrans(fname='SPHACC',nfreq=25,nmu=1,fmin=1e8,fmax=1e15,ename='SYNCHTHAV',nvals=1,spin=0.,mbh=1.,standard=1,nn=[10000,1,100],gridvals=[0.,400.,0.,0.],uout=.0025,oname='sphacc_abs.out')
            xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
        if save==0:
            i = pickle.load(open('test_grtrans_sphacc_intensity.p','rb'))
            if pgrtrans==0:
                terr = np.sum(np.abs(xlist[-1].ivals[:,0,14]-i))/np.sum(np.abs(i))
            else:
                terr = np.sum(np.abs(xlist[-1].ivals[0,:,14]-i))/np.sum(np.abs(i))
            print 'terr: ',terr
            if terr < (10*tol): passed+=1
            else: failed.append('sphacc intensity')
            max_passed+=1
            i = pickle.load(open('test_grtrans_sphacc_spectrum.p','rb'))
            terr = np.sum(np.abs(xlist[-1].spec-i))/np.sum(np.abs(i))
            print 'terr: ',terr
            if terr < (10*tol): passed+=1
            else: failed.append('sphacc spectrum')
            max_passed+=1
        else:
            pickle.dump(xlist[-1].ivals[:,0,14],open('test_grtrans_sphacc_intensity.p','wb'))
            pickle.dump(xlist[-1].spec,open('test_grtrans_sphacc_spectrum.p','wb'))

# ffjet
    xlist.append(gr.grtrans())
    if pgrtrans==0:
        xlist[-1].write_grtrans_inputs('inputs.in',fname='FFJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,400],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],ntscl=2.,nrscl=70.)
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='FFJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,400],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],ntscl=2.,nrscl=70.)
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_ffjet.p','rb'))
        if pgrtrans==0:
            terr = np.sum(np.abs(xlist[-1].ivals-i))/np.sum(np.abs(i))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals.transpose([1,0,2])-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('ffjet')
        max_passed+=1
    else:
        pickle.dump(xlist[-1].ivals,open('test_grtrans_ffjet.p','wb'))
# ffjet with delo integrator
    x2=gr.grtrans()
    if pgrtrans==0:
        x2.write_grtrans_inputs('inputs.in',fname='FFJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,1600],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],iname='delo',ntscl=2.,nrscl=70.)
        x2.run_grtrans()
        x2.read_grtrans_output()
    else:
        x2.run_pgrtrans(fname='FFJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,1600],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],iname='delo',ntscl=2.,nrscl=70.)
        x2.calc_spec_pgrtrans((np.shape(x2.ivals))[2])
    terr=10.
    terr = np.max(np.abs(x2.spec - xlist[-1].spec)/xlist[-1].spec)
    print 'terr: ',terr
    if terr < 0.05: passed += 1
    else: failed.append('delo')
    max_passed+=1
# ffjet with formal rad trans solution from Degl'Innocenti (1985):
    x3=gr.grtrans()
    if pgrtrans==0:
        x3.write_grtrans_inputs('inputs.in',fname='FFJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,1600],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],iname='formal',ntscl=2.,nrscl=70.)
        x3.run_grtrans()
        x3.read_grtrans_output()
    else:
        x3.run_pgrtrans(fname='FFJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,1600],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],iname='formal',ntscl=2.,nrscl=70.)
        x3.calc_spec_pgrtrans((np.shape(x3.ivals))[2])
    terr=10.
    terr = np.max(np.abs(x3.spec - xlist[-1].spec)/xlist[-1].spec)
    print 'terr: ',terr
    if terr < 0.05: passed += 1
    else: failed.append('formal')
    max_passed+=1
    x2=0; x3=0
# thindisk
    xlist.append(gr.grtrans())
    xlist[-1].write_grtrans_inputs('inputs.in',fname='THINDISK',nfreq=25,nmu=1,fmin=2.41e16,fmax=6.31e18,ename='BBPOL',nvals=4,spin=0.9,standard=2,nn=[100,100,1],uout=0.01,mbh=10, mumin=.26,mumax=.26,gridvals=[-21,21,-21,21])
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='THINDISK',nfreq=25,nmu=1,fmin=2.41e16,fmax=6.31e18,ename='BBPOL',nvals=4,spin=0.9,standard=2,nn=[100,100,1],uout=0.01,mbh=10, mumin=.26,mumax=.26,gridvals=[-21,21,-21,21])
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_thindisk.p','rb'))
        if pgrtrans==0:
            terr = np.sum(np.abs(xlist[-1].ivals-i))/np.sum(np.abs(i))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals.transpose([1,0,2])-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('thindisk')
        max_passed+=1
    else:
        pickle.dump(xlist[-1].ivals,open('test_grtrans_thindisk.p','wb'))
# total I w/, w/o pol
    xlist.append(gr.grtrans())
    xlist[-1].write_grtrans_inputs('inputs.in',fname='FFJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='SYNCHPL',nvals=1,spin=0.998,standard=1,nn=[100,100,400],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],ntscl=2.,nrscl=70.)
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='FFJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='SYNCHPL',nvals=1,spin=0.998,standard=1,nn=[100,100,400],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],ntscl=2.,nrscl=70.)
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_ffjet.p','rb'))
        if pgrtrans == 0:
            terr = np.sum(np.abs(xlist[-1].ivals[:,0,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals[0,:,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('polunpol')
        max_passed+=1

# harm
    xlist.append(gr.grtrans())
    xlist[-1].write_grtrans_inputs('inputs.in',fname='HARM',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=0.9375,standard=1,nn=[150,150,400],uout=0.04,mbh=4e6, mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,mumin=.6428,mumax=.6428,gridvals=[-13,13,-13,13],hhfile='dump040',hdfile='dump',hindf=40,hnt=1,muval=1./4.,gmin=1.)
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='HARM',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=0.9375,standard=1,nn=[150,150,400],uout=0.04,mbh=4e6, mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,mumin=.6428,mumax=.6428,gridvals=[-13.,13.,-13.,13.],fhfile='dump040',fdfile='dump',findf=40,fnt=1,muval=1./4.,gmin=1.)
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_harm.p','rb'))
        xlist[-1].ivals = np.where(xlist[-1].ivals==xlist[-1].ivals,xlist[-1].ivals,np.zeros(np.shape(xlist[-1].ivals)))
        i = np.where(i==i,i,np.zeros(np.shape(i)))
        if pgrtrans==0:
            terr = np.sum(np.abs(xlist[-1].ivals[:,0,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals[0,:,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('harm')
        max_passed+=1
    else:
        pickle.dump(xlist[-1].ivals,open('test_grtrans_harm.p','wb'))

    xlist.append(gr.grtrans())
    xlist[-1].write_grtrans_inputs('inputs.in',fname='POWERLAW',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHTH',nvals=4,spin=0.,standard=1,nn=[200,200,1600],uout=0.00005,mbh=4e6, mumin=0.5,mumax=0.5,nrotype=1,gridvals=[1200.,4000.,0.,2.*np.pi],iname='lsoda',srin=3200.,srout=3300.,ntscl=5e11,sthin=-0.02,sthout=0.02,rcut=4000.,snscl=1e5,phi0=-0.5,sphiin=0.,gmin=1.)
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='POWERLAW',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHTH',nvals=4,spin=0.,standard=1,nn=[200,200,1600],uout=0.00005,mbh=4e6, mumin=0.5,mumax=0.5,nrotype=1,gridvals=[1200.,4000.,0.,2.*np.pi],iname='lsoda',srin=3200.,srout=3300.,ntscl=5e11,sthin=-0.02,sthout=0.02,rcut=4000.,snscl=1e5,phi0=-0.5,sphiin=0.,gmin=1.)
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_toroidalfield.p','rb'))
        xlist[-1].ivals = np.where(xlist[-1].ivals==xlist[-1].ivals,xlist[-1].ivals,np.zeros(np.shape(xlist[-1].ivals)))
        i = np.where(i==i,i,np.zeros(np.shape(i)))
#        if pgrtrans==0:
#            terr = np.sum(np.abs(xlist[-1].ivals[:,0,0]-i[0,:,0]))/np.sum(abs(i[0,:,0]))
#        else:
#            terr = np.sum(np.abs(xlist[-1].ivals[0,:,0]-i[0,:,0]))/np.sum(abs(i[0,:,0]))
        if pgrtrans==0:
            terr = np.sum(np.abs(xlist[-1].ivals-i))/np.sum(np.abs(i))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals.transpose([1,0,2])-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < (2*tol): passed+=1
        else: failed.append('toroidal')
        max_passed+=1
    else:
        pickle.dump(xlist[-1].ivals,open('test_toroidalfield.p','wb'))

    print 'tests total: ', max_passed
    print 'tests passed: ', passed
    print 'tests failed: ',failed        

    return passed, max_passed, failed, xlist
