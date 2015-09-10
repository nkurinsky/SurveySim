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
    quit()

for fname in sys.argv[1:]:
    print fname
    newdir=fname.split('.fits')[0]
    os.mkdir(newdir)

    clf()
    figure()
    output=OutputFile(fname)
    output.info()
    output.counts.plot(1,xmin=5,xmax=400)

    bins=numpy.array([23.8,37.5,58.9,85.9,166.2,374.1])
    counts=numpy.array([2.0e8,6.4e7,1.2e7,3.1e6,2.1e5,1.7e4])*((bins*1e-3)**2.5)
    scatter(bins,counts,marker='v',label="Oliver et. al. 2010",color='red')

    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 250 $\mu$m")
    legend(loc='lower right')
    ylim(1,1e5)
    savefig(newdir+'/counts_250')
    clf()
    output.counts.plot(2,xmin=5,xmax=400)

    bins=numpy.array([23.8,37.5,58.9,85.9,166.2,374.1])
    counts=numpy.array([1.1e8,3.5e7,5.3e6,1.1e6,6.2e4,4.7e3])*((bins*1e-3)**2.5)
    scatter(bins,counts,marker='v',label="Oliver et. al. 2010",color='red')

    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 350 $\mu$m")
    legend(loc='lower right')
    ylim(1,1e5)
    savefig(newdir+'/counts_350')

    clf()
    output.counts.plot(3,xmin=5,xmax=400)

    bins=numpy.array([23.8,37.5,58.9,85.9,166.2,374.1])
    counts=numpy.array([3.6e7,1.1e7,1.6e6,2.3e5,1.3e4,1.3e3])*((bins*1e-3)**2.5)
    scatter(bins,counts,marker='v',label="Oliver et. al. 2010",color='red')
    
    ylabel('dN/dS $\mathrm{Jy}^{-1.5}$')
    title("SPIRE 500 $\mu$m")
    legend(loc='lower right')
    ylim(1,1e5)
    savefig(newdir+'/counts_500')

    clf()
    output.saveImages(newdir+'/images')
    if(output.fit()):
        output.MCMC.saveFit(newdir+'/fit')
        output.MCMC.saveFit(newdir+'/lffit',mode='lf')
        output.MCMC.saveFit(newdir+'/sedfit',mode='sed')
    print ""
