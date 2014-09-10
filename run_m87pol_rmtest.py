import grtrans_batch as gr
import numpy as np

mu=np.arange(10)*.1

#for i in range(len(mu[1:2])):
x=gr.grtrans()
x.compile_grtrans()
x.write_grtrans_inputs('inputs.in',fname='"MB09"',nfreq=10,fmin=5e10,fmax=5e13,ename='"POLSYNCHTH"',nvals=4,spin=0.9375,standard=1,nn="150,150,400",uout=0.04,mbh=6.3e9,mdotmin=1e20,mdotmax=1e20,nmdot=1,mumin=0.906,mumax=0.906,nmu=1,gridvals="-13,13,-13,13",muval=1./6.,nt=1,dt=20.)
x.run_production('m87pol_mbq_rmtest_spec')
