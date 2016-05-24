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

from SurveySim.OutputFile import *
from SurveySim.ModelFile import keyPrint

#Prelims to make prettier colors
palette=sns.color_palette()

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def plotImage(image,interpolation='nearest',labelCbar=True,annotateXY=None,cmax=105):
    print image.name
    tmpimg=np.flipud(image.img)
    cmap=cm.Greys
    clim=[0,cmax]
    if(image.name == 'Residual'):
        cmap=cm.bwr
        clim=[-cmax,cmax]
    if(image.name == 'Model Diagnostic'):
        cmap=cm.Blues
    if(image.name=='Observation Diagnostic'):
        cmap=cm.Reds

    plt.imshow(tmpimg, interpolation=interpolation, cmap=cmap,extent=image.extent, aspect='auto')
    if(image.name == 'Residual'):
        gpts=(tmpimg > 0)
        print 'St dev of residual: ',np.std(tmpimg[gpts])
        std_str=np.str(np.std(tmpimg[gpts]))
        if(annotateXY != None):
            plt.annotate('st.dev='+std_str[:4],xy=(annotateXY[0]+.15,annotateXY[1]-0.2),fontsize=15)
    plt.clim(clim)
    cbar=plt.colorbar()
    if(labelCbar):
        cbar.set_label("Normalized Density")

def plotImages(obs,xrange=None,yrange=None,axis1_name=None,axis2_name=None,annotateXY=None,cmax=105):
    col=0
    for key in obs.images.keys():
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
        if(key == 'Residual'):
            obs.norm = MidpointNormalize(midpoint=0)
        if(key == 'Observation Diagnostic'):
            obs.norm = MidpointNormalize(midpoint=-150)

        plotImage(obs.images[key],labelCbar=False,annotateXY=annotateXY,cmax=cmax)
        if(col > 1):
            plt.ylabel("")
        if(xrange != None):
            plt.xlim(xrange[0],xrange[1])
        if(yrange != None):
            plt.ylim(yrange[0],yrange[1])
        plt.xlabel(axis1_name)
        plt.ylabel(axis2_name)
        if(annotateXY != None):
            ax.annotate(lab,xy=annotateXY,fontsize=20)
        else:
            plt.title(lab,fontsize=20)
    plt.tight_layout(True)

def showImages(obs,block=True,xrange=None,yrange=None,axis1_name=None,axis2_name=None,annotateXY=None,cmax=105):
    plt.figure(figsize=(12,4))
    plotImages(obs,xrange=xrange,yrange=yrange,axis1_name=axis1_name,axis2_name=axis2_name,annotateXY=annotateXY,cmax=cmax)
    plt.show(block=block)

toPlot='mips'
if(toPlot == 'mips'):
    ofile='/Users/kurinsky/Desktop/F_mips_output.fits'
    output=OutputFile(ofile)
    axis1_name=r'log($S_{24}$/mJy)'
    axis2_name=r'log($S_{250}/S_{24}$)'
    showImages(output,
               xrange=[-1.3,1.0],
               yrange=[0.5,3],
               axis1_name=axis1_name,
               axis2_name=axis2_name,
               annotateXY=(-0.7,3*0.85))

elif(toPlot == 'spire'):
    ofile='Final_output_2/E_spire_output.fits'
    output=OutputFile(ofile)
    axis1_name=r'log($S_{250}$/mJy)'
    axis2_name=r'log($S_{350}/S_{250}$)'
    showImages(output,
               xrange=[0.9,2.8],
               yrange=[-1.0,1.0],
               axis1_name=axis1_name,
               axis2_name=axis2_name,
               annotateXY=(1.35,0.65))

