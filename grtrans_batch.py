from __future__ import print_function
import os
import namelist as nm
import numpy as np
#import matplotlib.pyplot as plt
from astropy.io import fits
# f2py grtrans module
from pgrtrans import pgrtrans
from time import time
#plt.ion()

pcG = 6.67259e-8
pcc2 = 8.98755179e20
msun = 1.99e33
cmperkpc=3.086e21

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
        self.uout=0.04
        self.uin=1.
        self.rcut=1.
        self.nrotype=2
        self.gridvals=[-25.,25.,-25.,25.]
        self.nn=[100,100,1]
        self.i1=-1
        self.i2=-1
        self.fname='THINDISK'
        self.dt=5.
        self.nt=1
        self.nload=1
        self.ename='BB'
        self.mbh=10.
        self.nmdot=1
        self.mdotmin=1.5e15
        self.mdotmax=1.5e15
        self.nfreq=1
        self.fmin=1.e17
        self.fmax=3.e19
        self.muval=1./4.
        self.sigcut=1.e10
        self.gmin=100
        self.gmax=1e5
        self.p1=3.5
        self.p2=3.5
        self.epotherargs=[1.,1.]
        self.epcoefindx=[1,1,1,1,1,1,1]
        self.jetalpha=0.02
        self.stype='const'
        self.use_geokerr='T'
        self.nvals=1
        self.iname='lsoda'
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
        self.jdfile='m87bl09rfp10xi5a998fluidvars.bin'
        self.ndfile='phatdiskm8st25.bin'
        self.ntscl=1.
        self.nrscl=6.
        self.snscl=3e7
        self.snnthscl=8e4
        self.snnthp=2.9
        self.sbeta=10.
        self.snp=0.
        self.snt=0.
        self.sbl06=0
        self.srin=-1.
        self.srout=1e8
        self.sthin=-10.
        self.sthout=10.
        self.sphiin=0.
        self.sphiout=1e4
        self.hrspot=1.5
        self.hr0spot=6.0
        self.hn0spot=1e4
        self.ofile='grtrans.out'
        self.outfile=''
        self.ifile='inputs.in'
        self.hhfile='dump040'
        self.hdfile='dump'
        self.hnt=1
        self.hindf=40
        self.tdfile='fieldline'
        self.tgfile='dump0000rr2.bin'
        self.tnt=1
        self.tnfiles=1
        self.tindf=5206
        self.tjonfix=1
        self.toff=0
        self.tsim='thickdiskrr2'
        self.tdindf=1
        self.tmagcrit=0
        self.mjonfix=1
        self.mnt=1
        self.mnfiles=1
        self.mgfile='mb09dipolegrid.bin'
        self.mdfile='/Users/Jason/analyses/mb09/dipole_rout1000/fieldline'
        self.mindf=1500
        self.msim='dipole'
        self.extra=0
        self.debug=0
# fluid arguments
        self.fdfile='fieldline'
        self.fgfile='dump0000rr2.bin'
        self.fhfile='dump040'
        self.fnt=1
        self.fnfiles=1
        self.findf=5206
        self.fjonfix=1
        self.foffset=0
        self.fsim='thickdiskrr2'
        self.fdindf=1
        self.fmagcrit=0
        self.fscalefac=1.
        self.fnw=500
        self.fwmin=1e-4
        self.fwmax=1e4
        self.fnfreq_tab=100
        self.ffmin=3.33e16
        self.ffmax=9e19
        self.frmax=1e4
        self.fnr=500
        self.fsigt=0.4
        self.ffcol=1.7
        self.frspot=1.5
        self.fr0spot=6.0
        self.fn0spot=1e4
        self.ftscl=1.
        self.frscl=6.
        self.fmdot=0.1

    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)
        if self.i1 < 0:
            self.i1 = 1
        if self.i2 < 0:
            self.i2 = self.nn[0]*self.nn[1]
        if self.i1 > self.nn[0]*self.nn[1]-1:
            print('fixing input i1 out of bounds: ',self.i1,self.nn[0]*self.nn[1])
            self.i1 = self.nn[0]*self.nn[1]-1
        if self.i2 > self.nn[0]*self.nn[1]:
            print('fixing input i2 out of bounds: ',self.i2,self.nn[0]*self.nn[1])
            self.i2 = self.nn[0]*self.nn[1]
        if self.ename == 'MAXCOMP':
            self.nweights = len(self.epotherargs)-2
        else:
            self.nweights = len(self.epotherargs)-1
        self.nep = len(self.epotherargs)
        self.delta = self.epotherargs[0]

    def write(self,fname):
        argst=['mdot','mbh']
        argsharm=['dfile','hfile','nt','indf'] #AC
        argskoral=['dfile','hfile','nt','indf']
        argsmb09=['gfile','dfile','nt','nfiles','indf','jonfix','sim']
        argsthick=['gfile','dfile','nt','nfiles','indf','jonfix','offset','sim','dindf','magcrit']
        argsp=['nw','wmin','wmax','nfreq_tab','fmin','fmax','rmax','nr','sigt','fcol']
        argsj=['dfile']
        argsn=['dfile','tscl','rscl']
        argsh=['rspot','r0spot','n0spot']
        names=['geodata','fluiddata','emisdata','general','harm','analytic'] #!AC ? 
        args=['standard','mumin','mumax','nmu','phi0','spin','uout','uin','rcut','nrotype','gridvals','nn','i1','i2','extra','debug']
        nargs=[len(args)]
        args1=['fname','dt','nt','nload','nmdot','mdotmin','mdotmax','sigcut']
        nargs.append(len(args1))
        args2=['ename','mbh','nfreq','fmin','fmax','muval','gmin','gmax','p1','p2','jetalpha','stype','delta','nweights','coefindx']
        nargs.append(len(args2))
        args3=['use_geokerr','nvals','iname','cflag']
        nargs.append(len(args3))
        args4=['fdfile','fgfile','fhfile','fnt','fnfiles','findf','fjonfix','foffset','fsim','fdindf','fmagcrit','fscalefac']
        nargs.append(len(args4))
        args5=['fnw','fwmin','fwmax','fnfreq_tab','ffmin','ffmax','frmax','fnr','fsigt','ffcol','frspot','fr0spot','fn0spot','ftscl','frscl','fmdot','fnscl','fnnthscl','fnnthp','fbeta','fbl06','fnp','ftp','frin','frout','fthin','fthout','fphiin','fphiout']
        nargs.append(len(args5))
        args=args
        args.extend(args1)
        args.extend(args2)
        args.extend(args3)
        args.extend(args4)
        args.extend(args5)
        nnstr = ''; gridstr = ''; cindxstr = ''
        for i in range(len(self.nn[:-1])):
            nnstr += str(self.nn[i])+','
        nnstr += str(self.nn[-1])
        for i in range(len(self.epcoefindx[:-1])):
            cindxstr += str(self.epcoefindx[i])+','
        cindxstr += str(self.epcoefindx[-1])
        for i in range(len(self.gridvals[:-1])):
            gridstr += str(self.gridvals[i])+','
        gridstr += str(self.gridvals[-1])
#        nnstr = str(self.nn[0])+','+str(self.nn[1])+','+str(self.nn[2])
#        gridstr = str(self.gridvals[0])+','+str(self.gridvals[1])+','+str(self.gridvals[
        vals=[self.standard,self.mumin,self.mumax,self.nmu,self.phi0,self.spin,self.uout,self.uin,self.rcut,self.nrotype,gridstr,nnstr,self.i1,self.i2,self.extra,self.debug]
        vals.extend(["'"+self.fname+"'",self.dt,self.nt,self.nload,self.nmdot,self.mdotmin,self.mdotmax,self.sigcut])
        vals.extend(["'"+self.ename+"'",self.mbh,self.nfreq,self.fmin,self.fmax,self.muval,self.gmin,self.gmax,self.p1,self.p2,self.jetalpha,"'"+self.stype+"'",self.delta,self.nweights,cindxstr])
        vals.extend([self.use_geokerr,self.nvals,"'"+self.iname+"'",self.cflag])
        print(self.fname)
# default values for fluid arguments that can change depending on model, work with thickdisk
        self.fdfile = self.tdfile; self.fgfile = self.tgfile; self.fhfile = self.hhfile
        self.fnt = self.tnt; self.fnfiles = self.tnfiles; self.findf = self.tindf
        self.fjonfix = self.tjonfix; self.foffset = self.toff; self.fsim = self.tsim
        self.fdindf = self.tdindf; self.fmagcrit = self.tmagcrit
        if self.fname=='MB09':
            namesmb09=['mb09']
            valsmb09=["'"+self.mgfile+"'","'"+self.mdfile+"'",self.mnt,self.mnfiles,self.mindf,self.mjonfix,self.msim]
            nargsmb09=[len(argsmb09)]
        # assign mb09 fluid arguments
            self.fdfile = self.mdfile
            self.fgfile = self.mgfile
            self.fnt = self.mnt
            self.findf = self.mindf
            self.fjonfix = self.mjonfix
            self.fsim = self.msim
            self.fnfiles = self.mnfiles

        if self.fname=='HARM':
            namesharm=['harm']
            valsharm=["'"+self.hdfile+"'","'"+self.hhfile+"'",self.hnt,self.hindf]
            nargsharm=[len(argsharm)]
        # assign harm fluid arguments
            self.fdfile = self.hdfile
            self.fhfile = self.hhfile
            self.fnt = self.hnt
            self.findf = self.hindf

        if self.fname=='KORAL':
            nameskoral=['koral']
            valskoral=["'"+self.hdfile+"'","'"+self.hhfile+"'",self.hnt,self.hindf]
            nargskoral=[len(argsharm)] #AC
            self.fdfile = self.hdfile
            self.fhfile = self.hhfile
            self.fnt = self.hnt
            self.findf = self.hindf

        if self.fname in ['KORAL3D','KORAL3D_DISK','KORAL3D_TOPJET','KORAL3D_BOTJET']:
            nameskoral=['koral']
            valskoral=["'"+self.hdfile+"'","'"+self.hhfile+"'",self.hnt,self.hindf]
            nargskoral=[len(argsharm)] #AC
            self.fdfile = self.hdfile
            self.fhfile = self.hhfile
            self.fnt = self.hnt
            self.findf = self.hindf

        if self.fname=='KORALNTH':
            nameskoral=['koral']
            valskoral=["'"+self.hdfile+"'","'"+self.hhfile+"'",self.hnt,self.hindf]
            nargskoral=[len(argsharm)] #AC
            self.fdfile = self.hdfile
            self.fhfile = self.hhfile
            self.fnt = self.hnt
            self.findf = self.hindf

        if self.fname=='THICKDISK':
            namesthick=['thickdisk']
            valsthick=["'"+self.tgfile+"'","'"+self.tdfile+"'",self.tnt,self.tnfiles,self.tindf,self.tjonfix,self.toff,"'"+self.tsim+"'",self.tdindf,self.tmagcrit]
            nargsthick=[len(argsthick)]
        if self.fname=='HARM3D':
#             harm fluid arguments
            self.fdfile = self.hdfile
            self.fhfile = self.hhfile
            self.fnt = self.hnt
            self.findf = self.hindf
            self.fgfile = self.hgfile
        if self.fname=='THINDISK':
            namest=['thindisk']
            valst=[self.tmdot,self.mbh]
            nargst=[len(argst)]
#            nm.write_namelist('thindisk.in',namest,argst,valst,nargst)
        if self.fname=='SHELL':
            namest=['shell']
            valst=[self.tmdot,self.mbh]
            nargst=[len(argst)]
        elif self.fname=='PHATDISK':
            namesp=['phatdisk']
            namest=['thindisk']
            valst=[self.tmdot,self.mbh]
            valsp=[self.pnw,self.pwmin,self.pwmax,self.pnfreq_tab,self.pfmin,self.pfmax,self.prmax,self.pnr,self.psigt,self.pfcol]
            nargst=[len(argst)]
            nargsp=[len(argsp)]
#            nm.write_namelist('thindisk.in',namest,argst,valst,nargst)
#            nm.write_namelist('phatdisk.in',namesp,argsp,valsp,nargsp)
        elif self.fname=='FFJET':
            namesj=['ffjet']
            valsj=["'"+self.jdfile+"'"]
            nargsj=[len(argsj)]
#            nm.write_namelist('ffjet.in',namesj,argsj,valsj,nargsj)
            self.fdfile = self.jdfile
        elif self.fname=='NUMDISK':
            namesn=['numdisk']
            valsn=["'"+self.ndfile+"'",self.ntscl,self.nrscl]
            nargsn=[len(argsn)]
#            nm.write_namelist('numdisk.in',namesn,argsn,valsn,nargsn)
        elif self.fname=='HOTSPOT':
            namesh=['hotspot']
            valsh=[self.hrspot,self.hr0spot,self.hn0spot]
            nargsh=[len(argsh)]
#            nm.write_namelist('hotspot.in',namesh,argsh,valsh,nargsh)
        elif self.fname=='SCHNITTMAN':
            namesh=['hotspot']
            valsh=[self.hrspot,self.hr0spot,self.hn0spot]
            nargsh=[len(argsh)]
#            nm.write_namelist('hotspot.in',namesh,argsh,valsh,nargsh)
#        else:
#            print('ERROR -- Unrecognized fluid name in grtrans_batch: ', self.fname)
#        print(args)
#        print(vals)
#        print(nargs)
        vals.extend(["'"+self.fdfile+"'","'"+self.fgfile+"'","'"+self.fhfile+"'",self.fnt,self.fnfiles,self.findf,self.fjonfix,self.foffset,"'"+self.fsim+"'",self.fdindf,self.fmagcrit,self.fscalefac])
        vals.extend([self.pnw,self.pwmin,self.pwmax,self.pnfreq_tab,self.pfmin,self.pfmax,self.prmax,self.pnr,self.psigt,self.pfcol,self.hrspot,self.hr0spot,self.hn0spot,self.ntscl,self.nrscl,self.tmdot,self.snscl,self.snnthscl,self.snnthp,self.sbeta,self.sbl06,self.snp,self.snt,self.srin,self.srout,self.sthin,self.sthout,self.sphiin,self.sphiout])
        nm.write_namelist(fname,names,args,vals,nargs)

class grtrans:
    """run & analyze fortran grtrans code"""
    def run_grtrans(self,**kwargs):
        if self.ofile=='grtrans.out':
            try:
                fh = open(self.ofile)
                fh.close()
                print('deleting grtrans.out...')
                os.system('rm -f grtrans.out')
            except IOError as e:
                pass
        else:
            try:
                #fh = open(self.ofile)
                #fh.close()
                #oname=self.ofile+str(time())
                #self.set_grtrans_input_file(self.ifile,oname)

                fh = open(self.ofile)
                fh.close()
                print('deleting %s...' % self.ofile)
                os.system('rm -f %s' %self.ofile)

            except IOError as e:
                pass
        os.system('./grtrans')

    def run_production(self,dir,**kwargs):
        self.ofile=dir+'/grtrans.out'
        self.set_grtrans_input_file(self.ifile,self.ofile)
        os.system('mkdir '+dir)
        os.system('cp *.in '+dir+'/.')
        self.run_grtrans()

    def compile_grtrans(self,**kwargs):
        os.system('make')

    def compile_pgrtrans(self,**kwargs):
        os.system('make pgrtrans')

    def compile_phat(self,**kwargs):
        os.system('./compilephat')

    def run_phat(self,**kwargs):
        os.system('./phat')
    
    def write_grtrans_inputs(self,ifile,oname='grtrans.out',**kwargs):
        self.inputs=grtrans_inputs(**kwargs)
        self.inputs.write(ifile)
        self.set_grtrans_input_file(ifile,oname)
        print(oname)

    def set_grtrans_input_file(self,fname,oname):
        nm.write_namelist('files.in',['files'],['ifile','ofile'],['"'+fname+'"','"'+oname+'"'],[2])
        self.ifile=fname
        self.ofile=oname

    def write_phat_inputs(self,a1,a2,ofile):
        nm.write_namelist('spins.in',['spins'],['a','athin','ofile'],[a1,a2,'"'+ofile+'"'],[3])

    def print_input_list(self):
        print('printing grtrans argument list')
        print(self.inputs.standard,self.inputs.mumin,self.inputs.mumax,self.inputs.nmu,self.inputs.phi0,self.inputs.spin,self.inputs.uout,self.inputs.uin,self.inputs.rcut,self.inputs.nrotype,self.inputs.gridvals,self.inputs.nn,self.inputs.i1,self.inputs.i2,self.inputs.fname,self.inputs.dt,self.inputs.nt,self.inputs.nload,self.inputs.nmdot,self.inputs.mdotmin,self.inputs.mdotmax,self.inputs.ename,self.inputs.mbh,self.inputs.nfreq,self.inputs.fmin,self.inputs.fmax,self.inputs.muval,self.inputs.gmin,self.inputs.gmax,self.inputs.p1,self.inputs.p2,self.inputs.jetalpha,self.inputs.stype,self.inputs.use_geokerr,self.inputs.nvals,self.inputs.iname,self.inputs.cflag,self.inputs.extra,self.inputs.debug,self.inputs.outfile,self.inputs.fdfile,self.inputs.fhfile,self.inputs.fgfile,self.inputs.fsim,self.inputs.fnt,self.inputs.findf,self.inputs.fnfiles,self.inputs.fjonfix,self.inputs.pnw,self.inputs.pnfreq_tab,self.inputs.pnr,self.inputs.foffset,self.inputs.fdindf,self.inputs.fmagcrit,self.inputs.hrspot,self.inputs.hr0spot,self.inputs.hn0spot,self.inputs.ntscl,self.inputs.nrscl,self.inputs.pwmin,self.inputs.pwmax,self.inputs.pfmin,self.inputs.pfmax,self.inputs.prmax,self.inputs.psigt,self.inputs.pfcol,self.inputs.tmdot,self.inputs.snscl,self.inputs.snnthscl,self.inputs.snnthp,self.inputs.sbeta,self.inputs.sbl06,self.inputs.snp,self.inputs.snt,self.inputs.srin,self.inputs.srout,self.inputs.sthin,self.inputs.sthout,self.inputs.sphiin,self.inputs.sphiout,self.inputs.fscalefac,self.inputs.sigcut,self.inputs.epcoefindx,self.inputs.epotherargs,self.inputs.nep)
        return

    def run_pgrtrans(self,**kwargs):
# do nothing if output data already exist
        if len(np.shape(pgrtrans.ivals)) > 0:
            print('ERROR in run_pgrtrans: data already exist!')
            return
# set arguments
        self.inputs=grtrans_inputs(**kwargs)
# call pgrtrans routine with arguments:
        self.print_input_list()
        pgrtrans.grtrans_main(self.inputs.standard,self.inputs.mumin,self.inputs.mumax,self.inputs.nmu,self.inputs.phi0,self.inputs.spin,self.inputs.uout,self.inputs.uin,self.inputs.rcut,self.inputs.nrotype,self.inputs.gridvals,self.inputs.nn,self.inputs.i1,self.inputs.i2,self.inputs.fname,self.inputs.dt,self.inputs.nt,self.inputs.nload,self.inputs.nmdot,self.inputs.mdotmin,self.inputs.mdotmax,self.inputs.ename,self.inputs.mbh,self.inputs.nfreq,self.inputs.fmin,self.inputs.fmax,self.inputs.muval,self.inputs.gmin,self.inputs.gmax,self.inputs.p1,self.inputs.p2,self.inputs.jetalpha,self.inputs.stype,self.inputs.use_geokerr,self.inputs.nvals,self.inputs.iname,self.inputs.cflag,self.inputs.extra,self.inputs.debug,self.inputs.outfile,self.inputs.fdfile,self.inputs.fhfile,self.inputs.fgfile,self.inputs.fsim,self.inputs.fnt,self.inputs.findf,self.inputs.fnfiles,self.inputs.fjonfix,self.inputs.pnw,self.inputs.pnfreq_tab,self.inputs.pnr,self.inputs.foffset,self.inputs.fdindf,self.inputs.fmagcrit,self.inputs.hrspot,self.inputs.hr0spot,self.inputs.hn0spot,self.inputs.ntscl,self.inputs.nrscl,self.inputs.pwmin,self.inputs.pwmax,self.inputs.pfmin,self.inputs.pfmax,self.inputs.prmax,self.inputs.psigt,self.inputs.pfcol,self.inputs.tmdot,self.inputs.snscl,self.inputs.snnthscl,self.inputs.snnthp,self.inputs.sbeta,self.inputs.sbl06,self.inputs.snp,self.inputs.snt,self.inputs.srin,self.inputs.srout,self.inputs.sthin,self.inputs.sthout,self.inputs.sphiin,self.inputs.sphiout,self.inputs.fscalefac,self.inputs.sigcut,self.inputs.epcoefindx,self.inputs.epotherargs,self.inputs.nep)
# read output
        self.ivals = pgrtrans.ivals.copy()
        self.ab = pgrtrans.ab.copy()
        self.nu = pgrtrans.freqs.copy()
        self.nx = int(self.inputs.nn[0])
        self.ny = int(self.inputs.nn[1])
        self.nvals = int(self.inputs.nvals)
        self.del_pgrtrans_data()

    def del_pgrtrans_data(self):
        if len(np.shape(pgrtrans.ivals))==0:
            print('ERROR in del_pgrtrans_data: no data to delete!')
            return
        pgrtrans.del_pgrtrans_data()

    def read_grtrans_output(self,bin=0):
        if bin==0:
            # fits read
            hdu=fits.open(self.ofile)
            n=len(hdu)
            print(len(hdu[0].data))
            ab=np.reshape(hdu[0].data,(int(len(hdu[0].data)/2),2))
# read image headers
            nu=np.empty(n-1); nx=np.empty(n-1); ny=np.empty(n-1)
            for i in range(n-1):
                nx[i]=hdu[i+1].header.get(6)
                ny[i]=hdu[i+1].header.get(7)
                nu[i]=hdu[i+1].header.get(8)

# assume each image has same nx, ny
            nx=int(nx[0]); ny=int(ny[0])
#            nx=int(np.sqrt(len(hdu[0].data)/2))
#            ny=int(nx)
            nvals=int(len(hdu[1].data)/nx/ny)
#            print('empty: ',nx,ny,nvals,n,np.shape(nx))
            ivals=np.empty((int(nx*ny),nvals,n-1))
            print(np.shape(ivals), np.shape(hdu[1].data))
# read images
            for i in range(n-1):
                ivals[:,:,i]=np.reshape(hdu[i+1].data,(int(nx*ny),int(nvals)))
                del hdu[i+1].data

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
        
        
        self.ab=ab
        self.ivals=ivals
        self.nu=nu
        self.nx=nx
        self.ny=ny
        self.nvals=nvals
        self.calc_spec(n-1)
#        self.calc_freqs(self.inputs.nfreq)

#!AC is this always right - no n?
    def calc_freqs(self,n):
        fmin=self.inputs.fmin
        fmax=self.inputs.fmax
   
        if(n>1):
            self.freqs=fmin*np.exp(np.arange(n)*np.log(fmax/fmin)/(n-1))
        else:
            self.freqs=np.array([fmin])
 
    def calc_spec(self,n):
#        print('test',np.shape(self.ab),np.shape(self.ivals))
#        if self.inputs.nrotype==1:
        if self.ny != 1:
            da=self.ab[self.nx,0]-self.ab[0,0]
            db=self.ab[1,1]-self.ab[0,1]
            spec=np.sum(self.ivals,0)*da*db
            if self.nvals>=4:
                self.lp=np.sqrt(spec[1,:]**2.+spec[2,:]**2.)/spec[0,:]
                self.cp=spec[3,:]/spec[0,:]
                self.lpf=np.sum(np.sqrt(self.ivals[:,1,:]**2.+self.ivals[:,2,:]**2.),0)*da*db/spec[0,:]
                self.cpf=np.sum(np.abs(self.ivals[:,3,:]),0)*da*db/spec[0,:]
        else:
            da=self.ab[1,0]-self.ab[0,0]
            db=0.
            spec=np.empty([n,int(self.inputs.nvals)])
            for i in range(n):
                for j in range(int(self.inputs.nvals)):
                    spec[i,j]=np.sum(self.ivals[:,j,i]*self.ab[:,0],0)*da*2.*3.14
    
        self.spec=spec
        self.da=da; self.db=db
#        self.convert_to_lum()


    def calc_spec_pgrtrans(self,n):
#       version for pgrtrans where array ordering is different
        if self.ny != 1:
            da=self.ab[0,self.nx]-self.ab[0,0]
            db=self.ab[1,1]-self.ab[1,0]
            spec=np.sum(self.ivals,1)*da*db
            if self.nvals>=4:
                self.lp=np.sqrt(spec[1,:]**2.+spec[2,:]**2.)/spec[0,:]
                self.cp=spec[3,:]/spec[0,:]
                self.lpf=np.sum(np.sqrt(self.ivals[1]**2.+self.ivals[2]**2.),0)*da*db/spec[0,:]
                self.cpf=np.sum(np.abs(self.ivals[3]),0)*da*db/spec[0,:]
        else:
            da=self.ab[0,1]-self.ab[0,0]
            db=0.
            spec=np.empty([n,int(self.inputs.nvals)])
            for i in range(n):
                for j in range(int(self.inputs.nvals)):
                    spec[i,j]=np.sum(self.ivals[j,:,i]*self.ab[0,:],0)*da*2.*np.arccos(-1.)
        self.spec = spec
#        self.convert_to_lum()
    
    def convert_to_lum(self):
        lbh = pcG*msun*self.inputs.mbh / pcc2 #!AC added msun
        fac= (4*np.pi*lbh**2)
        da = self.ab[self.nx,0]-self.ab[0,0]
        db = self.ab[1,1]-self.ab[0,1]
        print(da*db)
        self.spec *= fac
        self.ivals *= fac*da*db

# convert spectrum and intensities to L_iso in cgs units
    def convert_to_Jy(self,D):
        lbh = pcG*msun*self.inputs.mbh/pcc2 #!AC added msun
        #lbh = lbh/cmperkpc
        fac = (lbh**2/D**2.)*1e23
        da = self.ab[self.nx,0]-self.ab[0,0]
        db = self.ab[1,1]-self.ab[0,1]
        self.ivals *= fac*da*db
        self.spec *= fac

# test code 6/15/2019
# calculate image first and second moments and use to get xy centroids and semi-major/minor axes
    def calc_centroid_size(self,pgrtrans=-1):
        self.calc_spec(self.nx)
        nimages=len(self.spec[0])
        M00=np.sum(self.ivals[:,0],0)
        M10=np.zeros(nimages); M01=np.zeros(nimages); M11=np.zeros(nimages)
        M20=np.zeros(nimages); M02=np.zeros(nimages)
        for k in range(nimages):
            w=self.ivals[:,0,k]
            M10[k]=np.sum(w*self.ab[:,0])
            M01[k]=np.sum(w*self.ab[:,1])
            M20[k]=np.sum(w*self.ab[:,0]**2.)
            M02[k]=np.sum(w*self.ab[:,1]**2.)
            M11[k]=np.sum(w*self.ab[:,0]*self.ab[:,1])
        xcen=M10/M00; ycen=M01/M00; mu20=M20/M00-xcen**2.; mu11=M11/M00-xcen*ycen
        mu02=M02/M00-ycen**2.
# now do semi-major/minor axis and orientation in terms of moments
        theta=1./2.*np.arctan(2.*mu11/(mu20-mu02))
        fac=np.sqrt(4.*mu11**2.+(mu20-mu02)**2.)
        amax=np.sqrt((mu20+mu02+fac)/2.)
        amin=np.sqrt((mu20+mu02-fac)/2.)
        self.xcen=xcen; self.ycen=ycen; self.theta=theta
        self.amin=amin; self.amax=amax

    def disp_grtrans_image(self,idex=0,stokes=0,logscale=False,vmin=False,vmax=False):
        f=plt.figure()
        if logscale:
            ivals=self.ivals
            if vmin:
                ivals[ivals==0.] = 10**vmin
            else:
                ivals[ivals==0.] = 1.e-100 * np.mean(ivals)
            mask=ivals<=0
            ivals[mask] = np.min(ivals[~mask])
            ivals = np.log10(ivals)
        else:
            ivals = self.ivals



        if self.ivals.ndim < 3:
            imgplot = plt.imshow(np.transpose(ivals[:,stokes].reshape((self.nx,self.ny))),origin='lower',cmap='inferno',interpolation='gaussian')
        else:
            imgplot = plt.imshow(np.transpose(ivals[:,stokes,idex].reshape((self.nx,self.ny))),origin='lower',cmap='inferno',interpolation='gaussian')
        plt.axis('off')
        if vmin:
            imgplot.set_clim(vmin=vmin,vmax=vmax)
        plt.colorbar(imgplot, fraction=0.046, pad=0.04)
        plt.show()

    def disp_pgrtrans_image(self,idex=0,stokes=0):
        imgplot = plt.imshow(np.transpose(self.ivals[stokes,:,idex].reshape((self.nx,self.ny))),origin='lower')
        plt.show()

    def get_pol_vectors(self,idex=0,pgrtrans=1,nsamp=8,trim=-1,pcut=0.05):
        if trim < 0:
            nx=self.nx
            ny=self.ny
        else:
            nx=trim
            ny=trim
        if (np.mod(nx,nsamp)) != 0:
            for i in range(nsamp):
                nsamptry = nsamp-(i+1)
                print('nsamp change: ',i,nsamp-(i+1),np.mod(nx,nsamptry))
                if np.mod(nx,nsamptry)==0:
                    nsamp = nsamptry
                    break
            print('warning: nx/nsamp not an integer, changing nsamp = ',nsamp)
        X = np.arange(nx/nsamp,dtype=int)*nsamp+nsamp/2
        Y = np.arange(ny/nsamp,dtype=int)*nsamp+nsamp/2
        U,V = np.meshgrid(X,Y)
        if pgrtrans==1:
            evpa = 0.5*np.arctan2(self.ivals[2,:,idex],self.ivals[1,:,idex])
            m = np.sqrt(self.ivals[1,:,idex]**2.+self.ivals[2,:,idex]**2.)
            img = self.ivals[0,:,idex]
        else:
            evpa = 0.5*np.arctan2(self.ivals[:,2,idex],self.ivals[:,1,idex])
            m = np.sqrt(self.ivals[:,1,idex]**2.+self.ivals[:,2,idex]**2.)
            img = self.ivals[:,0,idex]
        scale=np.max(m)*10.
        mx = (np.transpose(np.resize(m * np.cos(evpa),(self.ny,self.nx))))
        my = (np.transpose(np.resize(m * np.sin(evpa),(self.ny,self.nx))))
        mask=np.zeros_like(mx)+1.; mask[np.resize(m/np.max(m),(self.ny,self.nx)) < pcut]=np.nan
        mx*=mask; my*=mask
        mx=mx[self.nx//2-nx//2:self.nx//2+nx//2,self.ny//2-ny//2:self.ny//2+ny//2]
        my=my[self.nx//2-nx//2:self.nx//2+nx//2,self.ny//2-ny//2:self.ny//2+ny//2]
        img=img.reshape((self.ny,self.nx)).transpose()
        img=img[self.nx//2-nx//2:self.nx//2+nx//2,self.ny//2-ny//2:self.ny//2+ny//2]
        print('img shape: ',np.shape(img),nx,self.nx)
        my=my[nsamp//2::nsamp,nsamp//2::nsamp]; mx=mx[nsamp//2::nsamp,nsamp//2::nsamp]
        return U,V,mx,my,img,scale


    def disp_pol_map(self,idex=0,pgrtrans=1,nsamp=8,sat=0.8,trim=-1,xlim=[-1.],ylim=[-1.],vmax=-1.,interp='bilinear',pcut=0.05):
        ###----------------------------------------------
#        i=img/np.max(img)/sat
        fig,ax = plt.subplots()
        left=0.16; bottom=0.16; right=0.95; top=0.96
        fig.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
        if np.sum(xlim)==-1.:
            ax.set_xlabel('alpha', fontsize=16)
            ax.set_ylabel('beta', fontsize=16)
            extent=[0,self.nx,0,self.nx]
        else:
            ax.set_xlabel('x (microarcseconds)',fontsize=14)
            ax.set_ylabel('y (microarcseconds)',fontsize=14)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            extent=[xlim[0],xlim[1],ylim[0],ylim[1]]
        leg = ax.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,frameon=False)
        self.pol_map(ax,idex=idex,pgrtrans=pgrtrans,nsamp=nsamp,sat=sat,trim=trim,extent=extent,vmax=vmax,interp=interp,pcut=pcut)

    def pol_map(self,ax,idex=0,pgrtrans=1,nsamp=8,sat=0.8,trim=-1,extent=[-1],vmax=-1.,interp='bilinear',pcut=0.05):
        if trim < 0:
            trim=self.nx
        if np.sum(extent)==-1:
            extent=[0,self.nx,0,self.nx]
        U,V,mx,my,img,scale = self.get_pol_vectors(idex=idex,pgrtrans=pgrtrans,nsamp=nsamp,trim=trim,pcut=pcut)
        if vmax < 0:
            vmax=sat*np.max(img)
#        img[img > sat*np.max(img)]=sat*np.max(img)
        ax.imshow(img,origin='lower',extent=extent,vmax=vmax,interpolation=interp)
        quiveropts = dict(color='white',headlength=0, pivot='middle', scale=scale,
                         width=6e-3, headwidth=1,headaxislength=0) # common options
        Uscl = (U-np.min(U))*1./(np.max(U)-np.min(U))*(extent[1]-extent[0])+extent[0]
        Vscl = (V-np.min(V))*1./(np.max(V)-np.min(V))*(extent[3]-extent[2])+extent[2]
        ax.quiver(Uscl,Vscl,mx,my,**quiveropts)

    def junhan_pol_diagnostics(self,ii=0,idex=0,pgrtrans=-1,nsamp=6,trim=128,fov=-1.,jyunit=60.):
# set field of view based on camera
# CHECK pgrtrans ordering implementing this for grtrans
        if fov==-1.:
            if pgrtrans==-1:
                fov=self.ab[-1,1]-self.ab[0,1]
            else:
                fov=self.ab[1,-1]-self.ab[1,0]
        if trim==-1:
            trim=self.nx
        plt.clf()
        fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
        if pgrtrans==-1:
#            self.calc_spec(self.nx)
            flux_tmp=[self.spec[0,idex]*jyunit,np.sqrt(self.spec[1,idex]**2.+self.spec[2,idex]**2.)/self.spec[0,idex]*100.,self.spec[3,idex]/self.spec[0,idex]*100.,0.5*np.arctan2(self.spec[2,idex],self.spec[1,idex])*180./np.pi]
            print('flux: ',flux_tmp)
        else:
            self.calc_spec_pgrtrans(self.nx)
            flux_tmp=[self.spec[0,idex]*jyunit,np.sqrt(self.spec[1,idex]**2.+self.spec[2,idex]**2.)/self.spec[0,idex]*100.,self.spec[3,idex]/self.spec[0,idex]*100.,0.5*np.arctan2(self.spec[2,idex],self.spec[1,idex])*180./np.pi]
        
        print(len(flux_tmp))
        plt.suptitle('Frame: #{0}, Frequency: {5:.0f} GHz\nFlux: {1:.1f} Jy,  '
                 'LP: {2:.1f}%,  '
                 'CP: {3:.1f}%, '
                 'EVPA: {4:.1f} deg'.
                     format(ii, flux_tmp[0], flux_tmp[1], flux_tmp[2], flux_tmp[3],self.nu[idex]/1e9))
        fov=fov*trim/self.nx*5.
        pol=['I','Q','U','V','P','LP','CP','EVPA']
        xlabel='x (microarcseconds)'; ylabel='y (microarcseconds)'
        if pgrtrans==-1:
            for k in range(4):
                ax = axes.flat[k]; ax.imshow(self.ivals[:,k,idex].reshape((self.nx,self.nx)).transpose()[self.nx/2-trim/2:self.nx/2+trim/2,self.nx/2-trim/2:self.nx/2+trim/2],origin='lower',extent=[-fov/2,fov/2,-fov/2,fov/2])
                ax.set_title(pol[k])
                ax.set_xlabel(xlabel)
                if k==0:
                    ax.set_ylabel(ylabel)
            ax = axes.flat[6]; ax.imshow((np.abs(self.ivals[:,3,idex])).reshape((self.nx,self.nx)).transpose()[self.nx/2-trim/2:self.nx/2+trim/2,self.nx/2-trim/2:self.nx/2+trim/2],origin='lower',extent=[-fov/2,fov/2,-fov/2,fov/2])
            ax.set_title(pol[6])
            ax.set_xlabel(xlabel)
            #ax.set_ylabel(ylabel)
            ax = axes.flat[5]; ax.imshow((np.sqrt(self.ivals[:,1,idex]**2.+self.ivals[:,2,idex]**2.)).reshape((self.nx,self.nx)).transpose()[self.nx/2-trim/2:self.nx/2+trim/2,self.nx/2-trim/2:self.nx/2+trim/2],origin='lower',extent=[-fov/2,fov/2,-fov/2,fov/2])
            ax.set_title(pol[5]); ax.set_xlabel(xlabel)#; ax.set_ylabel(ylabel)
            ax = axes.flat[4]; self.pol_map(ax,idex=idex,pgrtrans=pgrtrans,nsamp=nsamp,trim=trim,extent=[-fov/2,fov/2,-fov/2,fov/2])
            ax.set_title(pol[4]); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax = axes.flat[7]; ax.imshow((np.arctan2(self.ivals[:,2,idex],self.ivals[:,1,idex])).reshape((self.nx,self.nx)).transpose(),extent=[-fov/2,fov/2,-fov/2,fov/2],origin='lower')
            ax.set_title(pol[7]); ax.set_xlabel(xlabel)#; ax.set_ylabel(ylabel)
        else:
            for k in range(4):
                ax = axes.flat[k]; ax.imshow(self.ivals[k,:,idex].reshape((self.nx,self.nx)).transpose()[self.nx/2-trim/2:self.nx/2+trim/2,self.nx/2-trim/2:self.nx/2+trim/2],origin='lower',extent=[-fov/2,fov/2,-fov/2,fov/2])
                ax.set_title(pol[k]); ax.set_xlabel(xlabel)
                if k==0:
                    ax.set_ylabel(ylabel)
            ax = axes.flat[6]; ax.imshow((np.abs(self.ivals[3,:,idex])).reshape((self.nx,self.nx)).transpose()[self.nx/2-trim/2:self.nx/2+trim/2,self.nx/2-trim/2:self.nx/2+trim/2],origin='lower',extent=[-fov/2,fov/2,-fov/2,fov/2])
            ax.set_title(pol[6]); ax.set_xlabel(xlabel)#; ax.set_ylabel(ylabel)
            ax = axes.flat[5]; ax.imshow((np.sqrt(self.ivals[1,:,idex]**2.+self.ivals[2,:,idex]**2.)).reshape((self.nx,self.nx)).transpose()[self.nx/2-trim/2:self.nx/2+trim/2,self.nx/2-trim/2:self.nx/2+trim/2],origin='lower',extent=[-fov/2,fov/2,-fov/2,fov/2])
            ax.set_title(pol[5]); ax.set_xlabel(xlabel)#; ax.set_ylabel(ylabel)
            ax = axes.flat[4]; self.pol_map(ax,idex=idex,pgrtrans=pgrtrans,nsamp=nsamp,trim=trim)
            ax.set_title(pol[4]); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax = axes.flat[7]; ax.imshow((np.arctan2(self.ivals[2,:,idex],self.ivals[1,:,idex])).reshape((self.nx,self.nx)).transpose(),extent=[-fov/2,fov/2,-fov/2,fov/2],origin='lower')
            ax.set_title(pol[7]); ax.set_xlabel(xlabel)#; ax.set_ylabel(ylabel)
