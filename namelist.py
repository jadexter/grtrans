def write_namelist(fname,names,args,vals,nargs):
    f=open(fname,'w')
    indx=0
    for i in range(len(names)):
# this is the loop over namelist names
        f.write('&'+names[i]+'\n')
        for j in range(nargs[i]):
            f.write(' '+args[indx]+'='+str(vals[indx])+',\n')
            indx+=1

        f.write('/\n')

    f.close()
 
