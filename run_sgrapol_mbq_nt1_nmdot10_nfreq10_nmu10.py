import grtrans_batch as gr
import numpy as np

mu=np.arange(10)*.1

for i in range(len(mu[1:2])):
    x=gr.grtrans()
    x.compile_grtrans()
    x.write_grtrans_inputs('inputs.in',fname='"MB09"',nfreq=10,fmin=5e10,fmax=1e12,ename='"POLSYNCHTH"',nvals=4,spin=0.9375,standard=1,nn="150,150,400",uout=0.04,mbh=4e6, mdotmin=0.4e15,mdotmax=2e15,nmdot=5,mumin=mu[1],mumax=mu[1],nmu=1,gridvals="-13,13,-13,13",muval=1./4.,nt=1,dt=20.)
    x.run_production('sgrapol_mbq_m10_f10_'+str(mu[1]))
