#from pgriter import *

import pgriter

#example of using Alwin Mao's Iterator for Jason Dexter's GRTRANS

#in pgriter.py make sure putdir is defined, or comment out the line "from local import putdir"
#make sure modules sys, os, grtrans_batch, pgrtrans, and numpy are available

#This example uses sys.argv (command line arguments)
#example call:
#python pgrface.py 0.707106781187 0.0 1.0 2.0 0.25



alinputs,filename = pgriter.makealinputs()

iteramdots,iterafluxs,summaryvalue = pgriter.runiterator(alinputs,filename)

alinputs.snscl = iteramdots[pgriter.n.argmin(pgriter.n.abs(pgriter.n.array(iterafluxs) - pgriter.mastervalue))]

summaryname = pgriter.putdir + 'pgrsummary.txt'
summarytext = filename + '\t' + summaryvalue + '\n'
with open(summaryname,'a') as tfile:
    tfile.write(summarytext)

pgriter.makesariafspectra(filename,alinputs,summaryvalue)
