import grtrans_batch as gr
import pickle
import numpy as np
import copy

def run_test_problems(save=0,pgrtrans=0):
    # run grtrans test problems
    tol=1e-3; failed=[]
    x=gr.grtrans()
    if pgrtrans==0:
        x.compile_grtrans()
    else:
        x.compile_pgrtrans()
    passed=0; max_passed=0
# sphacc
#    x.write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=15,nmu=1,fmin=2.41e10,fmax=6.31e14,ename='POLSYNCHTH',nvals=4,spin=0.,mbh=10.,standard=1,uout=.003,nn="100,100,100")
# New tests of 1d intensity profile & full spectrum 12/14/2012
    x.write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=25,nmu=1,fmin=1e8,fmax=1e15,ename='SYNCHTHAV',nvals=1,spin=0.,mbh=1.,standard=1,nn="10000,1,100",gridvals="0.,400.,0.,0.",uout=.0025,oname='sphacc_abs.out')
    if pgrtrans==0:
        x.run_grtrans()
        x.read_grtrans_output()
    else:
        x.run_pgrtrans(fname='SPHACC',nfreq=25,nmu=1,fmin=1e8,fmax=1e15,ename='SYNCHTHAV',nvals=1,spin=0.,mbh=1.,standard=1,nn=[10000,1,100],gridvals=[0.,400.,0.,0.],uout=.0025,oname='sphacc_abs.out')
        x.calc_spec_pgrtrans((np.shape(x.ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_sphacc_intensity.p','rb'))
        if pgrtrans==0:
            terr = np.sum(np.abs(x.ivals[:,0,14]-i))/np.sum(np.abs(i))
        else:
            terr = np.sum(np.abs(x.ivals[0,:,14]-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('sphacc intensity')
        max_passed+=1
        i = pickle.load(open('test_grtrans_sphacc_spectrum.p','rb'))
        terr = np.sum(np.abs(x.spec-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('sphacc spectrum')
        max_passed+=1
    else:
        pickle.dump(x.ivals[:,0,14],open('test_grtrans_sphacc_intensity.p','wb'))
        pickle.dump(x.spec,open('test_grtrans_sphacc_spectrum.p','wb'))

    xlist = [x]
# toyjet
    xlist.append(gr.grtrans())
    if pgrtrans==0:
        xlist[-1].write_grtrans_inputs('inputs.in',fname='TOYJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn="100,100,400",uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals="-40,20,-20,40")
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='TOYJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,400],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40])
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_toyjet.p','rb'))
        if pgrtrans==0:
            terr = np.sum(np.abs(xlist[-1].ivals-i))/np.sum(np.abs(i))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals.transpose([1,0,2])-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('toyjet')
        max_passed+=1
    else:
        pickle.dump(xlist[-1].ivals,open('test_grtrans_toyjet.p','wb'))
# toyjet with delo integrator
    x2=gr.grtrans()
    if pgrtrans==0:
        x2.write_grtrans_inputs('inputs.in',fname='TOYJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn="100,100,1600",uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals="-40,20,-20,40",iname='delo')
        x2.run_grtrans()
        x2.read_grtrans_output()
    else:
        x2.run_pgrtrans(fname='TOYJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn=[100,100,1600],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40],iname='delo')
        x2.calc_spec_pgrtrans((np.shape(x2.ivals))[2])
    terr=10.
    terr = np.max(np.abs(x2.spec - xlist[-1].spec)/xlist[-1].spec)
    print 'terr: ',terr
    if terr < 0.05: passed += 1
    else: failed.append('delo')
    max_passed+=1
# thindisk
    xlist.append(gr.grtrans())
    xlist[-1].write_grtrans_inputs('inputs.in',fname='THINDISK',nfreq=25,nmu=1,fmin=2.41e16,fmax=6.31e18,ename='BBPOL',nvals=4,spin=0.9,standard=2,nn="100,100,1",uout=0.01,mbh=10, mumin=.26,mumax=.26,gridvals="-21,21,-21,21")
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
    xlist[-1].write_grtrans_inputs('inputs.in',fname='TOYJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='SYNCHPL',nvals=1,spin=0.998,standard=1,nn="100,100,400",uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals="-40,20,-20,40")
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='TOYJET',fdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='SYNCHPL',nvals=1,spin=0.998,standard=1,nn=[100,100,400],uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals=[-40,20,-20,40])
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_toyjet.p','rb'))
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
    xlist[-1].write_grtrans_inputs('inputs.in',fname='HARM',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=0.9375,standard=1,nn="150,150,400",uout=0.04,mbh=4e6, mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,mumin=.6428,mumax=.6428,gridvals="-13,13,-13,13",hhfile='dump040',hdfile='dump',hindf=40,hnt=1,muval=1./4.)
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='HARM',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=0.9375,standard=1,nn=[150,150,400],uout=0.04,mbh=4e6, mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,mumin=.6428,mumax=.6428,gridvals=[-13.,13.,-13.,13.],fhfile='dump040',fdfile='dump',findf=40,fnt=1,muval=1./4.)
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

# thickdisk
    xlist.append(gr.grtrans())
    xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=-0.9375,standard=1,nn="150,150,400",uout=0.04,mbh=4e6, mdotmin=0.5e13,mdotmax=0.5e13,nmdot=1,mumin=.906,mumax=.906,gridvals="-15,15,-15,15",tgfile='dump0000rr2.bin',tdfile='fieldline',tindf=5206,tnt=1,muval=1./41.,toff=0)
    if pgrtrans==0:
        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()
    else:
        xlist[-1].run_pgrtrans(fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=-0.9375,standard=1,nn=[150,150,400],uout=0.04,mbh=4e6, mdotmin=0.5e13,mdotmax=0.5e13,nmdot=1,mumin=.906,mumax=.906,gridvals=[-15,15,-15,15],fgfile='dump0000rr2.bin',fdfile='fieldline',findf=5206,fnt=1,muval=1./41.,foffset=0)
        xlist[-1].calc_spec_pgrtrans((np.shape(xlist[-1].ivals))[2])
    if save==0:
        i = pickle.load(open('test_grtrans_thick.p','rb'))
        if pgrtrans==0:
            terr = np.sum(np.abs(xlist[-1].ivals[:,0,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        else:
            terr = np.sum(np.abs(xlist[-1].ivals[0,:,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('thickdisk')
        max_passed+=1
    else:
        pickle.dump(xlist[-1].ivals,open('test_grtrans_thick.p','wb'))

    print 'tests total: ', max_passed
    print 'tests passed: ', passed
    print 'tests failed: ',failed

    return passed, max_passed, failed, xlist
