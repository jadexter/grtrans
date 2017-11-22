import os
from run_grtrans_test_problems import run_test_problems
from unit_tests import run_unit_tests

passed, max_passed, failed = run_test_problems(save=0)
nfailed, ufailed = run_unit_tests()
if passed < max_passed or nfailed > 0: print 'ERROR -- grtrans tests failed!'
else: 
#    os.chdir('..')
#    os.system('cvs -d :ext:jdexter@grad16.phys.washington.edu:/phys/users/jdexter/cvs commit grtrans')
    os.system('git commit -a')
    

#try:
#    with open('tests_failed.p','rb') as f: print 'ERROR -- grtrans tests failed!'
#except IOError as e:
#    os.chdir('..')
#    os.system('cvs -d :ext:jdexter@grad16.phys.washington.edu:/phys/users/jdexter/cvs commit grtrans')
