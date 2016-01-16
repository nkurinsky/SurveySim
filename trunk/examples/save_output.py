#!/usr/bin/env python

import os
import sys
import datetime
import time

sys.path.append('/usr/local/surveysim/python')

from OutputFile import *
from pylab import *
import numpy

rc('text', usetex=True)

if (len(sys.argv) < 2):
    print "Please provide model file name"
    quit()

for fname in sys.argv[1:]:
    print fname
    newdir=fname.split('.fits')[0]
    os.mkdir(newdir)

    clf()
    figure()
    output=OutputFile(fname)
    output.MCMC.info()
    output.saveImages(newdir+'/images')
    if(output.fit()):
        output.MCMC.saveFit(newdir+'/fit')
        output.MCMC.saveFit(newdir+'/lffit',mode='lf')
        output.MCMC.saveFit(newdir+'/sedfit',mode='sed_short')
    print ""
