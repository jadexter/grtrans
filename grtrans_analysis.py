import os
from locate import locate

def grtrans_thickdisk_analysis(simname,t0,tf):
    bpath="/lustre/medusa/jmckinne/data*/"
# make list of files
    for i in locate("fieldline*.bin",root=bpath+simname):
        for i in locate.locate("fieldline*.bin",root="/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/thickdisk7/"):
