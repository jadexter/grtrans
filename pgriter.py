#Alwin Mao's Iterator for grtrans as of June 23, 2015

#Added Secant method to the iterator, no longer producing 3 per iteration
#Factor of 4 speedup from before
#When using, look at the power law dependence of the flux on mdot or whatever, and use appropriately. 

#note iterator can be turned off by finding 33 + tries and replacing with -1 + tries 


import sys
import os 
from local import putdir
#this line is if you have a file local.py with directory information. only putdir is needed, to tell the iterator where to save files.  
import grtrans_batch as gr
from pgrtrans import pgrtrans
import numpy as n


x = gr.grtrans()
mecc = 5.92989e9
mastervalue = 3.4

#the value you are trying to match

#edit this function to define what constitutes a good enough fit
def isgoodenough(errors):
    three = len(n.where(errors < .01)[0]) > 3
    two = len(n.where(errors < .001)[0]) > 2
    one = len(n.where(errors < .0001)[0]) > 1
    return (three or two or one)

#least squares polynomial fit, avoids importing scipy...
def jgpoly(x,y,degree):
    array = n.zeros([len(x),degree])
    for no in n.arange(degree):
        array[:,no] = x**no
    matrix = n.mat(array)
    pseudo = (((matrix.T).dot(matrix)).I).dot((matrix.T).dot((n.mat(y)).T))
    wooodo = n.array(pseudo.reshape(degree))[0][::-1]
    return wooodo

#checks the new guess (newcenter) given the array of tested parameters (marr) and fluxes (farr), and the argument of the best run (best)
 
def checkguess(newcenter,marr,farr,best):
    bestflux = farr[best]
    bestmdot = marr[best]
    quality = True
    if True:
        if (bestflux < mastervalue) & (newcenter < bestmdot ):
            #Guess is shit. Flux is too small, we need bigger mdot
            quality = False
        elif (bestflux > mastervalue) & (newcenter > bestmdot):
            #Guess is shit. Flux is too big, we need a small mdot
            quality = False
        else:
            quality = True

#Sign check. Should use 0.0 but usually mdot shouldn't become so small.

    if newcenter < 1.0:
        quality = False

#make sure the new mass is actually new
    if (marr == newcenter).sum() > 0:
        quality = False

#don't jump more than 3 orders of magnitude at a time in case something funky happens
#but if you need to, allow jumping by 1 order of magnitude
#edit these if you want more/less flexibility
    if n.abs(n.log10(newcenter/bestmdot)) > 3:
        quality = False
#but please jump more than a part in a thousand...
    if n.abs(newcenter/bestmdot - 1.0) < 5e-3:
        quality = False

    return quality

#given the list of mdots and resulting fluxes, try many strategies to determine the next mdot to try.
#general strategy: try a safe method, check that it works and is a new guess, 
#if it doesn't lead to a better fit, checkguess prevents repeats, which forces finder() to try a riskier guess 

def finder(mdots,fluxs,tries=0,wander=2.0):
    sarr = n.argsort(mdots)
    marr = n.array(mdots)[sarr]
    farr = n.array(fluxs)[sarr]
    darr = ((farr - mastervalue)/mastervalue)**2.
    poss = n.where(farr > 0)
    best = n.where(darr == n.min(((farr[poss] - mastervalue)/mastervalue)**2.))[0][0] #only search through positive fluxs for the best
#    best = n.argmin(darr)
    quality = False

    #add data when there is not enough data or when run is just beginning

    if quality == False:
        if (len(darr) < (3 + tries)):
            newcenter = marr[best] * wander
            quality = checkguess(newcenter,marr,farr,best)
    if quality == False:
        if (len(darr) < (3 + tries)):
            newcenter = marr[best] / wander
            quality = checkguess(newcenter,marr,farr,best)


    if quality == False:
#if the mastervalue is in between 2 other guesses, use reliable midpoint methods
        upperbounds = (n.where(farr[poss] > mastervalue)[0])
        if len(upperbounds) > 0:
            lowerbounds = (n.where(farr[poss] < mastervalue)[0])
            if len(lowerbounds) > 0:
                #newcenter = n.sqrt(marr[upperbounds[0]]*marr[lowerbounds[-1]])
                newcenter = 0.5*(marr[poss][upperbounds[0]] + marr[poss][lowerbounds[-1]])
                quality = checkguess(newcenter,marr,farr,best)
                print('Testing Newcenter' + str(newcenter) + str(quality))
                if quality == False:
                    newcenter = n.sqrt(marr[upperbounds[0]]*marr[lowerbounds[-1]])
                    quality = checkguess(newcenter,marr,farr,best)
            else:
                quality = False
        else:
            quality = False

    #forced jump every 4 unless already done
#    if (marr == newcenter).sum() > 0:
    if quality == False:
        if (len(darr)%6 == 1):
            newcenter = marr[best] * wander
            quality = ((marr == newcenter).sum() == 0)
    if quality == False:
        if (len(darr)%6 == 2):
            newcenter = marr[best] / wander
            quality = ((marr == newcenter).sum() == 0)

    if quality == False:
        #secant method
        if len(darr[poss]) > 1:
            sbestdarr = n.sort(darr[poss])[1]
            sbest = n.where(darr == sbestdarr)[0][0]
            newcenter = marr[best] - (farr[best] - mastervalue) * (marr[best] - marr[sbest])/(farr[best] - farr[sbest])
            quality = checkguess(newcenter,marr,farr,best)

    #relies on I propto M**1/2.

    
    if quality == False:
        newcenter = marr[best] * (mastervalue/farr[best])**2.
        quality = checkguess(newcenter,marr,farr,best)

#clauses work with the order-of-magnitude limit imposed by checkguess
        if quality == False:
            newcenter = marr[best] * 10.
            quality = checkguess(newcenter,marr,farr,best)
            if quality == False:
                newcenter = marr[best] * 0.1
                quality = checkguess(newcenter,marr,farr,best)

    #relies on I propto M**2.

    if quality == False:
        newcenter = marr[best] / n.sqrt(farr[best]/mastervalue)
        quality = checkguess(newcenter,marr,farr,best)
        if quality == False:
            newcenter = marr[best] * 10.
            quality = checkguess(newcenter,marr,farr,best)
            if quality == False:
                newcenter = marr[best] * 0.1
                quality = checkguess(newcenter,marr,farr,best)

    #general power law

    if quality == False:
        if len(darr) > 1:
            sbest = n.argsort(darr)[1]
            newcenter = marr[best] * (mastervalue/farr[best])**(n.log(marr[best]/marr[sbest])/n.log(farr[best]/farr[sbest]))
            quality = checkguess(newcenter,marr,farr,best)
            if quality == False:
                newcenter = marr[best] * 10.
                quality = checkguess(newcenter,marr,farr,best)
                if quality == False:
                    newcenter = marr[best] * 0.1
                    quality = checkguess(newcenter,marr,farr,best)
        


    if quality == False:
        if len(marr[poss]) > 2:
            fit = jgpoly(n.array(marr[poss]),n.array(darr[poss]),3)
            newcenter = -.5 * fit[1]/fit[0]
            quality = checkguess(newcenter,marr,farr,best)
    if quality == False:
        if len(marr[poss]) > 2:
            a,b,c = jgpoly(marr[poss],farr[poss],3)
            newcenter = (-b + n.sqrt(b**2. - 4.0 * a * (c-mastervalue)))/(2.0*a)
            quality = checkguess(newcenter,marr,farr,best)

    #desperately try to find good data when there is not enough good data
    if quality == False:
        if len(poss) < 3:
            newcenter = marr[best] * 2.0
            quality = checkguess(newcenter,marr,farr,best)
    if quality == False:
        if len(poss) < 3:
            newcenter = marr[best] * 0.5
            quality = checkguess(newcenter,marr,farr,best)

    if quality == False:
        #Give up
        print('I give up. Try a more reasonable initial guess.')
        return newcenter,newcenter,True

    newwidth = n.sqrt(1.0+n.min(n.abs((n.array(mdots) - newcenter)/newcenter)))

    if len(marr) > 33 + tries:
        #This has gone on long enough
        print('I give up. Try a more reasonable guess.')
        return newwidth,newcenter,True
#old comments regarding newwidth
#zooms in, makes width smaller, if center is inside. 
#zooms as needed if center is outside
    goodenough = isgoodenough(darr)
#number of points within 10%
    return newwidth,newcenter,goodenough

#startiter

#inputs for SARIAF and MAXJUTT
class inputs:
    def __init__(self):
        self.standard=1
        self.fname='SARIAF'
        self.ename='MAXJUTT'
        self.nmu=1
        self.mumin=0.7#0.996194698092
        self.mumax=0.7#0.996194698092
        self.nfreq=1
        self.fmin=2.3e11
        self.fmax=2.3e11
        self.phi0=-0.5
        self.mbh=4.0e6
        self.spin=0.9375
        self.uout=0.005
        self.uin=1.0
        self.nvals=1
        self.gridvals=[-60.,60.,-60.,60.]
        self.nn=[150,150,400]
        self.muval=0.25
        self.stype='const'
        self.ntscl=10*mecc
        self.snscl=3e7
        self.snnthscl=8e4
        self.snnthp=2.9
        self.sbeta=10.
        self.sbl06=0
        self.delta=2.0
        self.coefs=[0.10714295,0.4459917,1.34219997,2.17450422,1.36515333,0.52892829]
        self.epotherargs = [self.delta] + self.coefs
        self.epcoefindx = [1,1,1,1,1,1,1]
def submit(alinputs):
    epotherargs = [float(alinputs.delta)] + alinputs.coefs
    x.run_pgrtrans(standard=alinputs.standard,fname=alinputs.fname,ename=alinputs.ename,nmu=alinputs.nmu,mumin=alinputs.mumin,mumax=alinputs.mumax,nfreq=alinputs.nfreq,fmin=alinputs.fmin,fmax=alinputs.fmax,phi0=alinputs.phi0,mbh=alinputs.mbh,spin=alinputs.spin,uout=alinputs.uout,uin=alinputs.uin,nvals=alinputs.nvals,gridvals=alinputs.gridvals,nn=alinputs.nn,muval=alinputs.muval,stype=alinputs.stype,ntscl=alinputs.ntscl,snscl=alinputs.snscl,snnthscl=alinputs.snnthscl,snnthp=alinputs.snnthp,sbeta=alinputs.sbeta,sbl06=alinputs.sbl06,epotherargs=alinputs.epotherargs,epcoefindx=alinputs.epcoefindx)
    x.calc_spec_pgrtrans(0)
    return (x.spec)*3.41e13

def makealinputs():
    alinputs=inputs()
    
    mu = float(sys.argv[1])
    spin = float(sys.argv[2])
    te = float(sys.argv[3])
    delta = float(sys.argv[4])
    muval = float(sys.argv[5])

#delta = 2.0
#te = 10.0
#spin = 0.0 
#mu = 0.7
#muval = 0.25

    alinputs.delta = delta
    alinputs.epotherargs[0] = delta
    alinputs.ntscl = te * mecc #theta_e * m_e c^2/ k_b in kelvin
    alinputs.spin = spin
    alinputs.mumin = mu
    alinputs.mumax = mu
    alinputs.muval = muval
    alinputs.delta = delta
    filename = putdir+'d'+ '{:0.3f}'.format(float(alinputs.delta))+'te'+'{:0.3f}'.format(te)+'a'+'{:0.3f}'.format(alinputs.spin)+'i'+'{:0.3f}'.format(alinputs.mumin)+'mvl'+'{:0.3f}'.format(alinputs.muval)+'sarjutt.txt'
    return alinputs,filename

#for SARIAF mdot is just n instead. Here we are. 
def runiterator(alinputs,filename):
    iteragoodenough = False
    iteramdots = []
    iterafluxs = []
    iterarange = 2.0
    iteracenter = alinputs.snscl #mdot is now the scale n number density
    param = 0 
#Load previous data
    iteratries = 0
    if os.path.isfile(filename):
        print('Loading ' + filename)
        data = n.loadtxt(filename)
        if len(data.shape) > 1:
            filtre = (data[:,2] == 2.3e11) & (data[:,3] > 0) & (data[:,1] > 0)
            iteramdots = list(data[filtre][:,1])
            iterafluxs = list(data[filtre][:,3])
            iteratries = len(iteramdots)
            iterarange,iteracenter,iteragoodenough = finder(iteramdots,iterafluxs,tries=iteratries)
    while iteragoodenough == False:
    #computation step
        param += 1
        alinputs.snscl = iteracenter 
        outflux = float(submit(alinputs)) 
    #submit runs it in pgrtrans
        print(outflux)
        outstring=str(param)+'\t'+str(alinputs.snscl)+'\t'+str(alinputs.fmin)+'\t'
        outstring+=str(outflux)+'\n'
        if param > -1:
            with open(filename,'a') as tfile:
                tfile.write(outstring)
#                pgrtrans.del_pgrtrans_data()
        iteramdots.append(alinputs.snscl)
        iterafluxs.append(outflux)
    #Iteration step
        iterarange,iteracenter,iteragoodenough = finder(iteramdots,iterafluxs,tries=iteratries)
    if isgoodenough((n.array(iterafluxs)/mastervalue - 1.0)**2.):
        summaryvalue = 'good'
    else:
        summaryvalue = 'bad'
    return iteramdots,iterafluxs,summaryvalue
#########################################


def makesariafspectra(filename,alinputs,summaryvalue):
    specfile = filename[:-4]+'.spec'

    if ((summaryvalue == 'good')):
        outspecnu = []
        outspecflux = []
        if os.path.isfile(specfile):
            oldspec = n.loadtxt(specfile)
            oldspec = oldspec.reshape(oldspec.size/2,2)
        else:
            oldspec = n.array([[0,0],[0,0]])

        flist = [1.5e10,2.7e10,4.868e10,8.77e10,1.58e11,2.0e11,6.0e11,2.0e12,3.0e12,3.0e13,2.5e11,1.5e12,4.0e12,1.6e13,9.0e13,2.3e11,1.0e10,3.1e9,1.0e9,3.0e14]


        for i in range(len(flist)):
            if (oldspec[:,0] == flist[i]).sum() == 0:
                if flist[i] <= 1.0e10:
                    camsize = 500
                elif flist[i] > 1.58e11:
                    camsize = 60
                else:
                    camsize = 100
                alinputs.fmin = flist[i]
                alinputs.fmax = flist[i]
                alinputs.nfreq = 1
                alinputs.gridvals=[-camsize,camsize,-camsize,camsize]
                print('Running ' + str(flist[i]) + '\t' + specfile)
                fluxs = submit(alinputs)[0]
                nus = x.nu
                for j in range(len(nus)):
                    specstring=str(nus[j])+'\t'+str(fluxs[j])+'\n'
                    with open(specfile,'a') as pfile:
                        pfile.write(specstring)

#alinputs,filename = makealinputs()

#iteramdots,iterafluxs,summaryvalue = runiterator(alinputs,filename)

#alinputs.snscl = iteramdots[n.argmin(n.abs(n.array(iterafluxs) - mastervalue))]

#summaryname = putdir + 'pgrsummary.txt'
#summarytext = filename + '\t' + summaryvalue + '\n'
#with open(summaryname,'a') as tfile:
#    tfile.write(summarytext)

#makesariafspectra(filename,alinputs,summaryvalue)




