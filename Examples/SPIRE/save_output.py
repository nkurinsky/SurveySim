#!/usr/bin/env python

import os
import sys
import datetime
import time
from OutputFile import *
from pylab import *
import numpy

rc('text', usetex=True)

if (len(sys.argv) < 2):
    print "Please provide model file name"
else:
    file=sys.argv[1]
    output=OutputFile(file)
    output.info()
    output.counts.plot(1,xmin=5,xmax=400)

    bins=numpy.array([23.8,37.5,58.9,85.9,166.2,374.1])
    counts=numpy.array([2.0e8,6.4e7,1.2e7,3.1e6,2.1e5,1.7e4])*((bins*1e-3)**2.5)
    scatter(bins,counts,marker='v',label="Oliver et. al. 2010",color='red')

    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 250 $\mu$m")
    legend(loc='lower right')
    ylim(1,1e5)
    savefig('counts_250')
    clf()
    output.counts.plot(2,xmin=5,xmax=400)

    bins=numpy.array([23.8,37.5,58.9,85.9,166.2,374.1])
    counts=numpy.array([1.1e8,3.5e7,5.3e6,1.1e6,6.2e4,4.7e3])*((bins*1e-3)**2.5)
    scatter(bins,counts,marker='v',label="Oliver et. al. 2010",color='red')

    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 350 $\mu$m")
    legend(loc='lower right')
    ylim(1,1e5)
    savefig('counts_350')

    clf()
    output.counts.plot(3,xmin=5,xmax=400)

    bins=numpy.array([23.8,37.5,58.9,85.9,166.2,374.1])
    counts=numpy.array([3.6e7,1.1e7,1.6e6,2.3e5,1.3e4,1.3e3])*((bins*1e-3)**2.5)
    scatter(bins,counts,marker='v',label="Oliver et. al. 2010",color='red')
    
    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 500 $\mu$m")
    legend(loc='lower right')
    ylim(1,1e5)
    savefig('counts_500')

    clf()
    output.saveImages('images')
    if(output.fit()):
        output.MCMC.saveFit('fit')
        output.MCMC.saveFit('lffit',mode='lf')
        output.MCMC.saveFit('sedfit',mode='sed')
    print ""
