#!/usr/bin/env python
# adding AzTEC 1.1 mm filter info to FILTER_LIST & allfilters.dat

# python imports
import os
import numpy as np
from astropy.io import ascii 
from subprocess import call
from scipy import constants

c = constants.c   # need this for frequency conversion 

# directory 
fdir='/Users/Jed/SurveySim/trunk/filters/'  # i ran this from inside python dir.

filterlist=fdir+"FILTER_LIST"
filterdata=fdir+"allfilters.dat"

dfilterlist=fdir+"FILTER_LIST_default"
dfilterdata=fdir+"default_allfilters.dat"

# open and read files
f=open(filterlist,'r')
flines=f.readlines()
nfilt=len(flines)
f.close()

print "number of filters: ", nfilt

newfiltid=nfilt

# open files and begin writing on new lines
f=open(filterlist,'a')
f.write('\n')                   
fd=open(filterdata,'a')
fd.write('\n')            ## REMEMBER TO UNCOMMENT THESE FOR FINAL RUN

# ---------------------# 
# Begin adding filter  #
# ---------------------#

newfiltername='AzTEC_1.1'

# read in data transmission data 
newfile='/Users/Jed/AzTEC_data/data/aztec_bandpass_05B.dat'
newdata=ascii.read(newfile)

# create wavelength and transmission arrays
freqs = np.array(newdata['col1']) # in GHz
trans = np.array(newdata['col2'])
lams = c / freqs # in angstroms (descending order)

# reverse arrays so lams is in ascending order
trans=trans[::-1]
lams=lams[::-1]
freqs=freqs[::-1]

newfiltid=newfiltid+1

to_ad=str(newfiltid)+"      "+newfiltername
f.write(to_ad)

to_ad='# '+newfiltername+' filter'
fd.write(to_ad)

i=0
for freqs in freqs:
    fd.write('\n')
    line='      '+str(lams[i])+' '+str(trans[i])
    fd.write(line)
    i=i+1


f.close()
fd.close()
