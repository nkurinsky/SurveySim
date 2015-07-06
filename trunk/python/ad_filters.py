#!/usr/bin/env python

#Python prelims
import os
import numpy as np
from subprocess import call

wd=os.getcwd()
fdir=wd+"/trunk/filters/"

filterlist=fdir+"FILTER_LIST"
filterdata=fdir+"allfilters.dat"

dfilterlist=fdir+"FILTER_LIST_default"
dfilterdata=fdir+"default_allfilters.dat"

#doesn't work, have to see how to get that sorted
#call["cp %s %s",(dfilterlist,filterlist)]
#call["cp %s %s",(dfilterdata,filterdata)]

f=open(filterlist,"r")
flines=f.readlines()
nfilt=len(flines)
f.close()

print 'Existing number of filters:', nfilt
newfiltid=nfilt
f=open(filterlist,"a")
f.write('\n') #ensure that start on a new line
fd=open(filterdata,"a")
fd.write('\n')

#========================================================
#start adding filters
#--------------------------------------------------------
newfiltername='ALMA band3'
freqs=np.arange(70,130,0.5) #GHz
lams=(2.99*1.e9)/freqs #in Angstroms as we use by default
trans=freqs*0.0
trans[np.logical_and(freqs > 84,freqs < 119)]=1

#these reverse the arrays so "lams" is in ascending order
trans=trans[::-1]
lams=lams[::-1]
freqs=freqs[::-1]

newfiltid=newfiltid+1

to_ad=str(newfiltid)+"     "+newfiltername
f.write(to_ad)

to_ad="# "+newfiltername+" filter"
fd.write(to_ad)

i=0
for freq in freqs:
    fd.write('\n')
    line="      "+str(lams[i])+" "+str(trans[i])
    fd.write(line)
    i=i+1

newfiltername='ALMA band6'
freqs=np.arange(200,300,0.5) #GHz
lams=(2.99*1.e9)/freqs #in Angstroms as we use by default
trans=freqs*0.0
trans[np.logical_and(freqs > 211,freqs < 275)]=1

#these reverse the arrays so "lams" is in ascending order
trans=trans[::-1]
lams=lams[::-1]
freqs=freqs[::-1]

newfiltid=newfiltid+1
f.write('\n')
to_ad=str(newfiltid)+"     "+newfiltername
f.write(to_ad)

fd.write('\n')
to_ad="# "+newfiltername+" filter"
fd.write(to_ad)
i=0
for freq in freqs:
    fd.write('\n')
    line="      "+str(lams[i])+" "+str(trans[i])
    fd.write(line)
    i=i+1

newfiltername='ALMA band7'
freqs=np.arange(250,400,0.5) #GHz
lams=(2.99*1.e9)/freqs #in Angstroms as we use by default
trans=freqs*0.0
trans[np.logical_and(freqs > 275,freqs < 370)]=1

#these reverse the arrays so "lams" is in ascending order
trans=trans[::-1]
lams=lams[::-1]
freqs=freqs[::-1]

newfiltid=newfiltid+1
f.write('\n')
to_ad=str(newfiltid)+"     "+newfiltername
f.write(to_ad)

to_ad="# "+newfiltername+" filter"
fd.write('\n')
fd.write(to_ad)
i=0
for freq in freqs:
    fd.write('\n')
    line="      "+str(lams[i])+" "+str(trans[i])
    fd.write(line)
    i=i+1

f.close()
fd.close()
