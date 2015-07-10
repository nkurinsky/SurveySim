#!/usr/bin/env python

import os
import sys
import datetime
import time
#sys.path.append("/usr/local/surveysim/python/")
from OutputFile import *
from pylab import *

rc('text', usetex=True)

if (len(sys.argv) < 2):
    print "Please provide model file name"
else:
    file=sys.argv[1]
    output=OutputFile(file)
    output.info()
    output.counts.plot(1)
    xlim(1,1e3)
    ylim(10,1e5)
    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 250 $\mu$m")
    savefig('counts_250')
    clf()
    output.counts.plot(2)
    xlim(1,1e3)
    ylim(10,1e5)
    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 350 $\mu$m")
    savefig('counts_350')
    clf()
    output.counts.plot(3)
    xlim(1,1e3)
    ylim(10,1e5)
    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 500 $\mu$m")
    savefig('counts_500')
    clf()
    output.saveImages('images')
    if(output.fit()):
        output.MCMC.saveFit('fit')
    print ""
