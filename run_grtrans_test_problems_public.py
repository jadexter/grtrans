import grtrans_batch as gr
import pickle
import numpy as np

def run_test_problems(save=0):
    # run grtrans test problems
    tol=1e-3; failed=[]
    x=gr.grtrans()
    x.compile_grtrans()
    passed=0; max_passed=0
# sphacc
    x=gr.grtrans()
    x.write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=15,nmu=1,fmin=2.41e10,fmax=6.31e14,ename='POLSYNCHTH',nvals=4,spin=0.,mbh=10.,standard=1,uout=.003,nn="100,100,100")
# New tests of 1d intensity profile & full spectrum 12/14/2012
    x.write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=25,nmu=1,fmin=1e8,fmax=1e15,ename='SYNCHTHAV',nvals=1,spin=0.,mbh=1.,standard=1,nn="10000,1,100",gridvals="0.,400.,0.,0.",uout=.0025,oname='sphacc_abs.out')
    x.run_grtrans()
    x.read_grtrans_output()
    if save==0:
        i = pickle.load(open('test_grtrans_sphacc_intensity.p','rb'))
        terr = np.sum(np.abs(x.ivals[:,0,14]-i))/np.sum(np.abs(i))
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
# ffjet
    x=gr.grtrans()
    x.write_grtrans_inputs('inputs.in',fname='TOYJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='POLSYNCHPL',nvals=4,spin=0.998,standard=1,nn="100,100,400",uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals="-40,20,-20,40")
    x.run_grtrans()
    x.read_grtrans_output()
    if save==0:
        i = pickle.load(open('test_grtrans_ffjet.p','rb'))
 #       x.ivals = np.where(x.ivals==x.ivals,x.ivals,np.zeros(np.shape(x.ivals)))
 #       i = np.where(i==i,i,np.zeros(np.shape(i)))
        terr = np.sum(np.abs(x.ivals-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('ffjet')
        max_passed+=1
    else:
        pickle.dump(x.ivals,open('test_grtrans_ffjet.p','wb'))
# thindisk
    x=gr.grtrans()
    x.write_grtrans_inputs('inputs.in',fname='THINDISK',nfreq=25,nmu=1,fmin=2.41e16,fmax=6.31e18,ename='BBPOL',nvals=4,spin=0.9,standard=2,nn="100,100,1",uout=0.01,mbh=10, mumin=.26,mumax=.26,gridvals="-21,21,-21,21")
    x.run_grtrans()
    x.read_grtrans_output()
    if save==0:
        i = pickle.load(open('test_grtrans_thindisk.p','rb'))
 #       x.ivals = np.where(x.ivals==x.ivals,x.ivals,np.zeros(np.shape(x.ivals)))
 #       i = np.where(i==i,i,np.zeros(np.shape(i)))
        terr = np.sum(np.abs(x.ivals-i))/np.sum(np.abs(i))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('thindisk')
        max_passed+=1
    else:
        pickle.dump(x.ivals,open('test_grtrans_thindisk.p','wb'))
# total I w/, w/o pol
    x=gr.grtrans()
    x.write_grtrans_inputs('inputs.in',fname='TOYJET',jdfile='m87bl09rfp10xi5a998fluidvars.bin',nfreq=1,nmu=1,fmin=3.45e11,fmax=3.45e11,ename='SYNCHPL',nvals=1,spin=0.998,standard=1,nn="100,100,400",uout=0.01,mbh=3.4e9, mumin=.906,mumax=.906,gridvals="-40,20,-20,40")
    x.run_grtrans()
    x.read_grtrans_output()
    if save==0:
        i = pickle.load(open('test_grtrans_ffjet.p','rb'))
#        x.ivals = np.where(x.ivals==x.ivals,x.ivals,np.zeros(np.shape(x.ivals)))
#        i = np.where(i==i,i,np.zeros(np.shape(i)))
        terr = np.sum(np.abs(x.ivals[:,0,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('polunpol')
        max_passed+=1

# harm
    x=gr.grtrans()
    x.write_grtrans_inputs('inputs.in',fname='HARM',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=0.9375,standard=1,nn="150,150,400",uout=0.04,mbh=4e6, mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,mumin=.6428,mumax=.6428,gridvals="-13,13,-13,13",hhfile='dump040',hdfile='dump',hindf=40,hnt=1,muval=1./4.)
    x.run_grtrans()
    x.read_grtrans_output()
    if save==0:
        i = pickle.load(open('test_grtrans_harm.p','rb'))
#        x.ivals = np.where(x.ivals==x.ivals,x.ivals,np.zeros(np.shape(x.ivals)))
#        i = np.where(i==i,i,np.zeros(np.shape(i)))
        terr = np.sum(np.abs(x.ivals[:,0,0]-i[:,0,0]))/np.sum(abs(i[:,0,0]))
        print 'terr: ',terr
        if terr < tol: passed+=1
        else: failed.append('harm')
        max_passed+=1
    else:
        pickle.dump(x.ivals,open('test_grtrans_harm.p','wb'))

    print 'tests total: ', max_passed
    print 'tests passed: ', passed
    print 'tests failed: ',failed
#    if passed < max_passed: pickle.dump(failed,open('tests_failed.p','wb'))

    return passed, max_passed, failed
#run()
