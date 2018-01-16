# compile and run grtrans unit tests
# based on 12/21/2012 notes

import os
import numpy as np
import grtrans_batch as gr

def run_unit_tests(grtrans_dir='',compile=0):

    failed = []
    nfailed = 0

    if compile==1:
# compile grtrans library
        os.system('make all')


    if grtrans_dir=='':
        grtrans_dir = os.getcwd()
# kerr tests: parallel transport of pol, consistency check between two pol methods
    angtol=1e-2; dottol=1e-6
    x=gr.grtrans()
# write appropriate inputs for comparison to geokerr
    x.write_grtrans_inputs('inputs.in',fname='THINDISK',nfreq=25,nmu=1,fmin=2.41e16,fmax=6.31e18,ename='BBPOL',nvals=4,spin=0.9,standard=2,nn=[100,100,1],uout=0.01,mbh=10, mumin=.26,mumax=.26,gridvals=[-21,21,-21,21])
    os.system('rm unit_test_kerr.out')
    os.system('gfortran test_kerr.f90 -L'+grtrans_dir+' -lgrtrans')
    os.system('./a.out')
    maxang,maxangdiff,kdotk,kdota,perpk = np.loadtxt('unit_test_kerr.out')
    if maxang > angtol or maxangdiff > dottol or kdotk > dottol or kdota > dottol or perpk > dottol:
        print 'Error in kerr unit test!'
        failed.append('kerr')
        nfailed+=1
    
# geodesics unit tests: k^\mu ?= dx^\mu / d\lambda
    

# fluid unit tests
    fluid_tests = ['hotspot','harm','ffjet']
    ubtol = [1e-4, 1e-2, 1e-1, 0.4]; utol = [1e-4, 1e-2, 1e-1, 0.4]
    for i in range(len(fluid_tests)):
        print 'i: ',i,range(len(fluid_tests)),fluid_tests[i]
        os.system('rm unit_test_'+fluid_tests[i]+'.out')
        os.system('cp '+fluid_tests[i]+'.in.dist '+fluid_tests[i]+'.in')
        os.system('gfortran test_'+fluid_tests[i]+'.f90 -L'+grtrans_dir+' -lgrtrans')
        os.system('./a.out')
        minn, maxnorm, maxub = np.loadtxt('unit_test_'+fluid_tests[i]+'.out')
        print 'vals: ',minn, maxnorm, maxub
        if maxnorm > utol[i] or maxub > ubtol[i] or minn < 0:
            print 'Error in '+fluid_tests[i]+' unit test!'
            failed.append(fluid_tests[i])
            nfailed+=1
# integration tests

    if nfailed==0: print 'Passed all unit tests!'
    else: print 'Failed tests: ',failed
    print 'Exiting unit tests'

    return nfailed,failed
