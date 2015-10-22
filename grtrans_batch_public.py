import os
import namelist as nm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab
import pyfits
from time import time
pylab.ion()

def flatten(l):
    if isinstance(l,list):
        return sum(map(flatten,l))
    else:
        return l

class grtrans_inputs:
    def init(self):
        self.cflag=1
        self.standard=2
        self.mumin=.1
        self.mumax=1.
        self.phi0=-0.5
        self.nmu=10
        self.spin=0.998
        self.uout=0.0001
        self.uin=1.
        self.rcut=1.
        self.nrotype=2
        self.gridvals='-25.,25.,-25.,25.'
        self.nn='100,100,1'
        self.fname='"THINDISK"'
        self.dt=5.
        self.nt=1
        self.nload=1
        self.ename='"BB"'
        self.mbh=10.
        self.nmdot=1
        self.mdotmin=1.5e15
        self.mdotmax=1.5e15
        self.nfreq=1
        self.fmin=1.e17
        self.fmax=3.e19
        self.muval=1./4.
        self.gmin=100
        self.gmax=1e5
        self.p1=3.5
        self.p2=3.5
        self.use_geokerr='T'
        self.nvals=1
        self.iname='"lsoda"'
        self.tmdot=0.1
        self.pnw=500
        self.pwmin=1.e-4
        self.pwmax=1.e4
        self.pnfreq_tab=100
        self.pfmin=self.fmin/3.
        self.pfmax=self.fmax*3.
        self.prmax=1.e4
        self.pnr=500
        self.psigt=0.4
        self.pfcol=1.7
        self.jdfile='"m87bl09rfp10xi5a998fluidvars.bin"'
        self.ndfile='"phatdiskm8st25.bin"'
        self.ntscl=30.
        self.nrscl=6.
        self.hrspot=1.5
        self.hr0spot=6.0
        self.hn0spot=1e4
        self.ofile='grtrans.out'
        self.ifile='inputs.in'
        self.hhfile='"dump040"'
        self.hdfile='"dump"'
        self.hnt=1
        self.hindf=40

    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    def write(self,fname):
        argst=['mdot','mbh']
        argsharm=['dfile','hfile','nt','indf']
        argsp=['nw','wmin','wmax','nfreq_tab','fmin','fmax','rmax','nr','sigt','fcol']
        argsj=['dfile']
        argsn=['dfile','tscl','rscl']
        argsh=['rspot','r0spot','n0spot']
        names=['geodata','fluiddata','emisdata','general']
        args=['standard','mumin','mumax','nmu','phi0','spin','uout','uin','rcut','nrotype','gridvals','nn']
        nargs=[len(args)]
        args1=['fname','dt','nt','nload','nmdot','mdotmin','mdotmax']
        nargs.append(len(args1))
        args2=['ename','mbh','nfreq','fmin','fmax','muval','gmin','gmax','p1','p2']
        nargs.append(len(args2))
        args3=['use_geokerr','nvals','iname','cflag']
        nargs.append(len(args3))
        args=args
        args.extend(args1)
        args.extend(args2)
        args.extend(args3)
        vals=[self.standard,self.mumin,self.mumax,self.nmu,self.phi0,self.spin,self.uout,self.uin,self.rcut,self.nrotype,self.gridvals,self.nn]
        vals.extend([self.fname,self.dt,self.nt,self.nload,self.nmdot,self.mdotmin,self.mdotmax])
        vals.extend([self.ename,self.mbh,self.nfreq,self.fmin,self.fmax,self.muval,self.gmin,self.gmax,self.p1,self.p2])
        vals.extend([self.use_geokerr,self.nvals,self.iname,self.cflag])
        print self.fname
        if self.fname=='"HARM"':
            namesharm=['harm']
            valsharm=[self.hdfile,self.hhfile,self.hnt,self.hindf]
            nargsharm=[len(argsharm)]
#            print valsharm
#            print nargsharm
#            print argsharm
            nm.write_namelist('harm.in',namesharm,argsharm,valsharm,nargsharm)
        if self.fname=='"THINDISK"':
            namest=['thindisk']
            valst=[self.tmdot,self.mbh]
            nargst=[len(argst)]
            nm.write_namelist('thindisk.in',namest,argst,valst,nargst)
        elif self.fname=='"PHATDISK"':
            namesp=['phatdisk']
            namest=['thindisk']
            valst=[self.tmdot,self.mbh]
            valsp=[self.pnw,self.pwmin,self.pwmax,self.pnfreq_tab,self.pfmin,self.pfmax,self.prmax,self.pnr,self.psigt,self.pfcol]
            nargst=[len(argst)]
            nargsp=[len(argsp)]
            nm.write_namelist('thindisk.in',namest,argst,valst,nargst)
            nm.write_namelist('phatdisk.in',namesp,argsp,valsp,nargsp)
        elif self.fname=='"FFJET"':
            namesj=['ffjet']
            valsj=[self.jdfile]
            nargsj=[len(argsj)]
            nm.write_namelist('ffjet.in',namesj,argsj,valsj,nargsj)
        elif self.fname=='"NUMDISK"':
            namesn=['numdisk']
            valsn=[self.ndfile,self.ntscl,self.nrscl]
            nargsn=[len(argsn)]
            nm.write_namelist('numdisk.in',namesn,argsn,valsn,nargsn)
        elif self.fname=='"HOTSPOT"':
            namesh=['hotspot']
            valsh=[self.hrspot,self.hr0spot,self.hn0spot]
            nargsh=[len(argsh)]
            nm.write_namelist('hotspot.in',namesh,argsh,valsh,nargsh)
        elif self.fname=='"SCHNITTMAN"':
            namesh=['hotspot']
            valsh=[self.hrspot,self.hr0spot,self.hn0spot]
            nargsh=[len(argsh)]
            nm.write_namelist('hotspot.in',namesh,argsh,valsh,nargsh)
#        print args
#        print vals
#        print nargs
        nm.write_namelist(fname,names,args,vals,nargs)

class grtrans:
    """run & analyze fortran grtrans code"""
    def run_grtrans(self,**kwargs):
        if self.ofile=='grtrans.out':
            try:
                fh = open(self.ofile)
                fh.close()
                print 'deleting grtrans.out...'
                os.system('rm -f grtrans.out')
            except IOError as e:
                pass
        else:
            try:
                fh = open(self.ofile)
                fh.close()
                oname=self.ofile+str(time())
                self.set_grtrans_input_file(self.ifile,oname)
            except IOError as e:
                pass
        os.system('./grtrans')

    def compile_grtrans(self,**kwargs):
        os.system('make')

    def compile_phat(self,**kwargs):
        os.system('./compilephat')

    def run_phat(self,**kwargs):
        os.system('./phat')
    
    def write_grtrans_inputs(self,iname,oname='grtrans.out',**kwargs):
        self.inputs=grtrans_inputs(**kwargs)
        self.inputs.write(iname)
        self.set_grtrans_input_file(iname,oname)
        print oname

    def set_grtrans_input_file(self,fname,oname):
        nm.write_namelist('files.in',['files'],['ifile','ofile'],['"'+fname+'"','"'+oname+'"'],[2])
        self.ifile=fname
        self.ofile=oname

    def write_phat_inputs(self,a1,a2,ofile):
        nm.write_namelist('spins.in',['spins'],['a','athin','ofile'],[a1,a2,'"'+ofile+'"'],[3])

    def read_grtrans_output(self,bin=0):
        if bin==0:
            # fits read
            hdu=pyfits.open(self.ofile)
            n=len(hdu)
            print len(hdu[0].data)
            ab=np.reshape(hdu[0].data,(len(hdu[0].data)/2,2))
# read image headers
            nu=np.empty(n-1); nx=np.empty(n-1); ny=np.empty(n-1)
            for i in range(n-1):
                nx[i]=hdu[i+1].header.get(6)
                ny[i]=hdu[i+1].header.get(7)
                nu[i]=hdu[i+1].header.get(8)

# assume each image has same nx, ny
            nx=nx[0]; ny=ny[0]
#            nx=int(np.sqrt(len(hdu[0].data)/2))
#            ny=int(nx)
            nvals=len(hdu[1].data)/nx/ny
#            print 'empty: ',nx,ny,nvals,n,np.shape(nx)
            ivals=np.empty((nx*ny,nvals,n-1))
            print np.shape(ivals), np.shape(hdu[1].data)
# read images
            for i in range(n-1):
                ivals[:,:,i]=np.reshape(hdu[i+1].data,(nx*ny,nvals))

        else:
            with open(self.ofile,'rb') as f:
                temp,nx,ny,nvals=np.fromfile(f,dtype='i4',count=4)
                temp=np.fromfile(f,dtype='i4',count=2)
                nkey=np.fromfile(f,dtype='i4',count=1)
                temp=np.fromfile(f,dtype='i4',count=2)
                keys=np.fromfile(f,dtype='f4',count=nkey)
                temp=np.fromfile(f,dtype='i4',count=2)
                field=np.fromfile(f,dtype='f4',count=2*nx*ny)
                ab=np.reshape(field,(nx*ny,2))
                temp=np.fromfile(f,dtype='i4',count=2)
                field=np.fromfile(f,dtype='f4',count=nvals*nx*ny)
                intens=np.reshape(field,(nx*ny,nvals))
                ivals=intens
                nu=keys[0]
                for x in f:
                    temp,nx,ny,nvals=np.fromfile(x,dtype='i4',count=4)
                    temp=np.fromfile(x,dtype='i4',count=2)
                    nkey=np.fromfile(x,dtype='i4',count=1)
                    temp=np.fromfile(x,dtype='i4',count=2)
                    keys=np.fromfile(x,dtype='f4',count=nkey)
                    temp=np.fromfile(x,dtype='i4',count=2)
                    field=np.fromfile(x,dtype='f4',count=2*nx*ny)
                    ab=np.reshape(field,(nx*ny,2))
                    temp=np.fromfile(x,dtype='i4',count=2)
                    field=np.fromfile(x,dtype='f4',count=nvals*nx*ny)
                    intens=np.reshape(field,(nx*ny,nvals))
                    ivals=np.append(ivals,intens)
                    nu=nu.append(nu,keys[0])
        
        print 'test',np.shape(ab),np.shape(ivals)            
        if ny != 1:
            da=ab[nx,0]-ab[0,0]
            db=ab[1,1]-ab[0,1]
            spec=np.sum(ivals,0)*da*db
        else:
            da=ab[1,0]-ab[0,0]
            db=0.
            spec=np.empty([n-1,int(nvals)])
            for i in range(n-1):
                for j in range(int(nvals)):
                    spec[i,j]=np.sum(ivals[:,j,i]*ab[:,0],0)*da*2.*3.14
    
        self.ab=ab
        self.ivals=ivals
        self.nu=nu
        self.spec=spec
        self.nx=nx
        self.ny=ny

    def disp_grtrans_image(self,idex):
        if self.ivals.ndim < 3:
            imgplot = plt.imshow(np.transpose(self.ivals[:,0].reshape((self.nx,self.ny))),origin='lower')
        else:
            imgplot = plt.imshow(np.transpose(self.ivals[:,0,idex].reshape((self.nx,self.ny))),origin='lower')
        plt.show()
