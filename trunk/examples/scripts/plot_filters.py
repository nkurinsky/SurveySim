#!/usr/bin/env python

#Python prelims
import os
import sys
from pylab import *
import numpy
from SurveySim.filters import *

#my_filters=['Herschel_SPIRE_250um', 'Herschel_SPIRE_350um', 'Herschel_SPIRE_500um']
my_filters=['Spitzer_IRAC_1','Spitzer_IRAC_2','Spitzer_IRAC_4']

ids=[]
for filt in my_filters:
    fid,fname=getFilterID(filt)
    print(fname)
    ids.append(fid)
print(ids)

clf()
fdata=fill_filters(ids)
for i in range(0,3):
    x=numpy.array(fdata[i],dtype=numpy.float)*1e-10
    plot(x,fdata[i+3],label=my_filters[i])
xlabel('Wavelength')
ylabel('Transmission')
ylim(0,1.3)
xscale('log')
legend()
show()
