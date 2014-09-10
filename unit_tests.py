# compile and run grtrans unit tests
# based on 12/21/2012 notes

import os
import numpy as np
import grtrans_batch as gr

def run_unit_tests():

    failed = []
    nfailed = 0
    grtrans_dir='/Users/Jason/code/grtrans'

# compile grtrans as library
    os.system('make lib')

# first is inputs unit test
#    geotol=1e-3; geokerrtol=1e-6
# run geokerr for test case
#    os.system('geokerr < geokerr_unit_test1.in > geokerr_unit_test1.out')
#    os.system('geokerr < geokerr_unit_test2.in > geokerr_unit_test2.out')
#    x=gr.grtrans()
# write appropriate inputs for comparison to geokerr
#    x.write_grtrans_inputs()
#    os.system('rm unit_test_geo.out')
#    os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_geo.f90')
#    os.system('./a.out')
#    maxgeodiff,maxgeokerrdiff1,maxgeokerrdiff2 = np.loadtxt('unit_test_geo.out')
#    if maxgeodiff > geotol or maxgeokerrdiff1 > geokerrtol or maxgeokerrdiff2 > geokerrtol:
#        print 'Error in geodesics unit test!'
#        failed.append('geodesics')
#        nfailed+=1

# kerr tests: parallel transport of pol, consistency check between two pol methods
    angtol=1e-2; dottol=1e-6
    x=gr.grtrans()
# write appropriate inputs for comparison to geokerr
    x.write_grtrans_inputs('inputs.in',fname='THINDISK',nfreq=25,nmu=1,fmin=2.41e16,fmax=6.31e18,ename='BBPOL',nvals=4,spin=0.9,standard=2,nn="100,100,1",uout=0.01,mbh=10, mumin=.26,mumax=.26,gridvals="-21,21,-21,21")
    os.system('rm unit_test_kerr.out')
    os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_kerr.f90')
    os.system('./a.out')
    maxang,maxangdiff,kdotk,kdota,perpk = np.loadtxt('unit_test_kerr.out')
    if maxang > angtol or maxangdiff > dottol or kdotk > dottol or kdota > dottol or perpk > dottol:
        print 'Error in kerr unit test!'
        failed.append('kerr')
        nfailed+=1
    
# geodesics unit tests: k^\mu ?= dx^\mu / d\lambda
    

# fluid unit tests
    fluid_tests = ['hotspot','harm','toyjet','thickdisk']
    ubtol = [1e-4, 1e-2, 1e-1, 0.4]; utol = [1e-4, 1e-2, 1e-1, 0.4]
    for i in range(len(fluid_tests)):
        print 'i: ',i,range(len(fluid_tests)),fluid_tests[i]
        os.system('rm unit_test_'+fluid_tests[i]+'.out')
        os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_'+fluid_tests[i]+'.f90')
        os.system('./a.out')
        minn, maxnorm, maxub = np.loadtxt('unit_test_'+fluid_tests[i]+'.out')
        print 'vals: ',minn, maxnorm, maxub
        if maxnorm > utol[i] or maxub > ubtol[i] or minn < 0:
            print 'Error in '+fluid_tests[i]+' unit test!'
            failed.append(fluid_tests[i])
            nfailed+=1
#os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_numdisk.f90')
#os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_toyjet.f90')

# emissivity tests
#os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_emis.f90')
#os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_interpemis.f90')
#os.system('gfortran -L'+grtrans_dir+' -lgrtrans test_interpemis_phatdisk.f90')

# integration tests

    if nfailed==0: print 'Passed all unit tests!'
    else: print 'Failed tests: ',failed
    print 'Exiting unit tests'

    return nfailed,failed
