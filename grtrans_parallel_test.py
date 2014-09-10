import grtrans_batch as gr
import os
from time import time

#def grtrans_parallel_test():
x = gr.grtrans()
x.compile_grtrans()
x.write_grtrans_inputs('inputs.in',fname='SPHACC',nfreq=25,nmu=1,fmin=1e8,fmax=1e15,ename='SYNCHTHAV',nvals=1,spin=0.,mbh=1.,standard=1,nn="10000,1,100",gridvals="0.,400.,0.,0.",uout=.0025,oname='sphacc_abs.out')
threads = ['1','2','4','8','16','32','64']; t = []
for i in range(len(threads)):
#    os.system('export OMP_NUM_THREADS='+threads[i])
    start = time()
#    x.run_grtrans()
    os.system('./run_grtrans_parallel.sh '+threads[i])
    end = time()
    t.append(end-start)

print 'threads, time: ', threads,t
