#!/usr/bin/env python

#Python prelims
import os
import sys

from filters import *
import matplotlib.pyplot as plt
import numpy as np

instruments=["AzTEC"] #"WISE", "JWST", "Herschel_SPIRE" etc.

print(filterDir())

nfilt=0
for instrument in instruments:
    ids,names=getFilterIDs(instrument)
    print instrument,names

nfilt=len(ids)

if(nfilt == 3):
    lam1,lam2,lam3,trans1,trans2,trans3=fill_filters(ids)
#convert to micron from default Angstrom
    lam1=np.array(lam1)/10000.
    lam2=np.array(lam2)/10000.
    lam3=np.array(lam3)/10000.
    plt.plot(lam1,trans1,lam2,trans2,lam3,trans3)
if(nfilt == 4):
    lam1,lam2,lam3,lam4,trans1,trans2,trans3,trans4,=fill_filters(ids)
    lam1=np.array(lam1)/10000.
    lam2=np.array(lam2)/10000.
    lam3=np.array(lam3)/10000.
    lam4=np.array(lam4)/10000.
    plt.plot(lam1,trans1,lam2,trans2,lam3,trans3,lam4,trans4)

plt.xtitle='wavelength [um]'
plt.ytitle='Transmission'
plt.show()
