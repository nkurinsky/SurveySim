#!/usr/bin/env python

import pyfits as fits
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import numpy as np
from numpy import std
from scipy import stats
from scipy.stats import norm
from scipy.interpolate import interp1d
import math
from math import ceil
import os
import Tkinter as tk


print 'COSMOS'
obsfile='/usr/local/surveysim/obs/L2-COSMOS_xID24_DR3.fits'
        #Single-band catalog
#obsfile="/usr/local/surveysim/obs/L5-Lockman-SWIRE_xID250_DR2.fits"
#obsfile="/usr/local/surveysim/obs/L2-COSMOS_xID250_DR2.fits"
f = fits.open(obsfile) 
tbdata = f[1].data
ra=tbdata['RA']
dec=tbdata['DEC']
f250=tbdata['F250'] #if using band-merged or SExtractor

f24=tbdata['F24'] 
        ##f250=tbdata['Flux'] #if using single-band Sussex extraction
#gpts = (f250 > 8.0)

gpts = (f24 > 60.0)

ra_gpts=ra[gpts]
dec_gpts=dec[gpts]

#print len(ra_gpts)
plt.plot(ra_gpts,dec_gpts,'.')
plt.show()

def areacoverage(ra,dec,dra,ddec):
    n_srcs=len(ra)
    rmin=np.min(ra)
    rmax=np.max(ra)
    dmin=np.min(dec)
    dmax=np.max(dec)

    nra=ceil((rmax-rmin)/dra)
    ndec=ceil((dmax-dmin)/ddec)
    #print nra*ndec
    tile=np.zeros((nra,ndec))
    #print np.sum(tile)
    for i in range(0,n_srcs):
        side1=ceil((ra[i]-rmin)/dra)
        side2=ceil((dec[i]-dmin)/ddec)
        if((side1 < nra) and (side2 < ndec)):
            tile[side1,side2]=1

    nhits=np.sum(tile)
    decmean=np.mean(dec)
    tilearea=dra*ddec*np.cos(decmean*(math.pi/180.0))
    #print tilearea
    area=nhits*tilearea
    return area

area=areacoverage(ra_gpts,dec_gpts,0.05,0.05)

print area

#based on Fig12 in Bethermin et al. 
testarea=1.4*1.4*np.cos(2.2*(math.pi/180.0))
print testarea
