#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from ModelFile import *
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv) < 2):
    print "Please provide model file name"
else:
    file=sys.argv[1]
    mod=ModelFile()
    mod.load(file)
    mod.info()
    #print mod.filters[0].fid
    ids=[mod.filters[0].fid,mod.filters[1].fid,mod.filters[2].fid]
    print ids
    lam1,lam2,lam3,trans1,trans2,trans3=fill_filters(ids)
#convert to micron from default Angstrom
    lam1=np.array(lam1)/10000.
    lam2=np.array(lam2)/10000.
    lam3=np.array(lam3)/10000.
    plt.plot(lam1,trans1,lam2,trans2,lam3,trans3)
    plt.xtitle='wavelength [um]'
    plt.ytitle='Transmission'
    plt.show()

