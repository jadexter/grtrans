import grtrans_batch as gr

x=gr.grtrans()
x.compile_grtrans()
x.write_grtrans_inputs('inputs.in',fname='MB09',nfreq=3,fmin=2.3e11,fmax=6.9e11,ename='POLSYNCHTH',nvals=4,spin=0.94,standard=1,nn="150,150,400",uout=0.04,mbh=4e6, mdotmin=0.8e15,mdotmax=0.8e15,nmdot=1,mumin=0.6428,mumax=0.6428,nmu=1,gridvals="-13,13,-13,13",muval=1./6.,nt=1,dt=20.)
x.run_production('sgrapol_mbqt5_psphere_nofar')
