import grtrans_batch as gr
import numpy as np

x=gr.grtrans()
x.compile_grtrans()
x.write_grtrans_inputs('inputs.in',fname='MB09',nfreq=10,fmin=5e10,fmax=1e12,ename='POLSYNCHTH',nvals=4,spin=0.92,standard=1,nn="150,150,400",uout=0.04,mbh=4e6, mdotmin=1e15,mdotmax=7e15,nmdot=5,mumin=0.0,mumax=0.9,nmu=10,gridvals="-13,13,-13,13",muval=1./4.,nt=1,dt=20.)
x.run_production('sgrapol_mbd_m10_f10_mdot5')
