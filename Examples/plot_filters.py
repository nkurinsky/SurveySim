#!/usr/bin/env python

#Python prelims
import os
import sys
from pylab import *

sys.path.append("../Python/")
from filters import *

my_filters=['Herschel_SPIRE_250um', 'Herschel_SPIRE_350um', 'Herschel_SPIRE_500um']

ids=[]
for filt in my_filters:
    fid,fname=getFilterID(filt)
    print(fname)
    ids.append(fid)
print(ids)

fdata=fill_filters(ids)
for i in range(0,3):
    plot(fdata[i],fdata[i+3])
xlabel('Wavelength')
ylabel('Transmission')
xscale('log')
show()
