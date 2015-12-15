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
from math import ceil
import os
import Tkinter as tk
from matplotlib.colors import Normalize

#from OutputFile import *

#Prelims to make prettier colors
import seaborn as sns
sns.set(style='white',font='serif',font_scale=1.5,palette='Set2')
#sns.despine()
palette=sns.color_palette()

plt.rc('text', usetex=True)

params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{sfmath}']}
plt.rcParams.update(params)

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def keyPrint(dict):
    for key,val in dict.iteritems():
        print '\t',key,'-\t',val

def binSize (dat): 
    return 3.49*std(dat)/pow(len(dat),0.333)

def nbins (dat): 
    return int(ceil((max(dat)-min(dat))/binSize(dat)))

class FitInfo:

    exclude=['SIMPLE',
             'H_MIN_X','H_MAX_X','BINSIZE_X',
             'H_MIN_Y','H_MAX_Y','BINSIZE_Y',
             'FITSTAT']

    def __init__(self,hdus,fitted):
        phdr=hdus[0].header
        self.initial=dict()
        for key in phdr.keys():
            if not (key in FitInfo.exclude):
                self.initial[key] = phdr[key]


class FitImage:
    
    colorMaps={'Residual':cm.bwr,
              'model':cm.Greys,
              'observation':cm.Greens}

    def __init__(self,hdus,extension):
        self.img=hdus[extension].data

        if(extension == 0):
            self.img=hdus[2].data-hdus[1].data
            #print len(self.img)
            #print np.nonzero(self.img)
            #tmp=np.nonzero(self.img)
            #tmp1=np.absolute(tmp[2])
            #print tmp1
            #print np.mean(tmp1),np.median(tmp1)

        self.phdu=hdus[0].header

        self.xbinsize=self.phdu['BINSIZE_X']
        self.ybinsize=self.phdu['BINSIZE_Y']

        #print self.xbinsize,self.ybinsize
        self.xmin=self.phdu['H_MIN_X']
        self.xmax=self.phdu['H_MAX_X']
        self.ymin=self.phdu['H_MIN_Y']
        self.ymax=self.phdu['H_MAX_Y']

        self.extent=[self.ymin,self.ymax,self.xmin,self.xmax]
        if(extension == 0):
            self.name='Residual'
        else:
            self.name=hdus[extension].header['EXTNAME']


    def __str__(self):
        return "Image: "+self.name

    def plot(self,interpolation='nearest',labelCbar=True):
        print self.name
        tmpimg=np.flipud(self.img)
        cmap=cm.Greys
        clim=[0,105]
        if(self.name == 'Residual'):
            cmap=cm.bwr
            clim=[-105,105]
        if(self.name == 'Model Diagnostic'):
            cmap=cm.Blues
        if(self.name=='Observation Diagnostic'):
            cmap=cm.Reds

        plt.imshow(tmpimg, interpolation=interpolation, cmap=cmap,extent=self.extent, aspect='auto')
        plt.clim(clim)
        cbar=plt.colorbar()
        if(labelCbar):
            cbar.set_label("Normalized Density")

    def show(self):
        self.plot()
        plt.show(block=True)

    def save(self,imgfile):
        plt.figure()
        self.plot()
        savefig(imgfile)

class OutputFile:

    def __init__(self,filename):
        self.filename=filename
        self.hdus=fits.open(filename)
        self.hdus.info
        self.images=dict()
        if(len(self.hdus) == 3):
            self.fitExt=False
            self.runtype='Simulation Only'
            self.simInfo=FitInfo(self.hdus,False)
        else:
            self.fitExt=True
            self.runtype='Parameter Fitting'
            self.simInfo=FitInfo(self.hdus,True)
            for i in range(0,3):
                self.images['temp'] = FitImage(self.hdus,i)
                newkey=self.images['temp'].name
                self.images[newkey] = self.images.pop('temp')


    def plotImages(self,xrange=None,yrange=None,axis1_name=None,axis2_name=None):
        col=0
        for key in self.images.keys():
            if(key == "Observation Diagnostic"):
                col=1
                lab='Observed'
            if(key == "Residual"):
                col=3
                lab='Residual'
            if(key == "Model Diagnostic"):
                col=2
                lab='Model'
            ax=plt.subplot(1,3,col)
            #plt.subplot(1,3,col)
            if(key == 'Residual'):
                self.norm = MidpointNormalize(midpoint=0)
            if(key == 'Observation Diagnostic'):
                self.norm = MidpointNormalize(midpoint=-150)

            self.images[key].plot(labelCbar=False)
            if(col > 1):
                plt.ylabel("")
            if(xrange != None):
                plt.xlim(xrange[0],xrange[1])
            if(yrange != None):
                plt.ylim(yrange[0],yrange[1])
            plt.xlabel(axis1_name)
            plt.ylabel(axis2_name)
            #MIPS
            #ax.annotate(lab,xy=(-0.7,yrange[1]*0.85),fontsize=20)
            #SPIRE
            ax.annotate(lab,xy=(1.35,0.65),fontsize=20)
        plt.tight_layout(True)

    def showImages(self,block=True,xrange=None,yrange=None,axis1_name=None,axis2_name=None):
        plt.figure(figsize=(12,4))
        self.plotImages(xrange=xrange,yrange=yrange,axis1_name=axis1_name,axis2_name=axis2_name)
        plt.show(block=block)

    def saveImages(self,filename,xrange=None,yrange=None):
        plt.figure(figsize=(12,4))
        self.plotImages(xrange=xrange,yrange=yrange)
        plt.savefig(filename)


    def show(self,block=True):
        my_dpi=116 
        maxsize=maxWindowSize()
        figsize=(maxsize[0]/float(my_dpi)-0.25,maxsize[1]/float(my_dpi)-0.5)
        plt.figure(figsize=figsize, dpi=my_dpi)
        self.plot()
        plt.show(block=block)


ofile='I_diffdiagn_output.fits'
mfile='I_diffdiagn_model.fits'

output=OutputFile(ofile)

#MIPS
#axis1_name=r'log($S_{24}$/mJy)'
#axis2_name=r'log($S_{250}/S_{24}$)'
#SPIRE
axis1_name=r'log($S_{250}$/mJy)'
axis2_name=r'log($S_{350}/S_{250}$)'


#MIPS
#output.showImages(xrange=[-1.3,1.0],yrange=[0.5,3],axis1_name=axis1_name,axis2_name=axis2_name)
#SPIRE
output.showImages(xrange=[0.9,2.8],yrange=[-1.0,1.0],axis1_name=axis1_name,axis2_name=axis2_name)

