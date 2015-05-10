import pyfits as fits
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy
from numpy import std
from math import ceil

from ModelFile import keyPrint

def bins (dat): return int(ceil(pow(len(dat),0.333)))
def binSize (dat): return 2.4*std(dat)/pow(len(dat),0.333)
def nbins (dat): return int(ceil((max(dat)-min(dat))/binSize(dat)))
def nbins2 (dat): return int(ceil((max(dat)-min(dat))/binSize(dat)/2))
def contour2d (x,y):
    H, xedges, yedges = numpy.histogram2d(y, x, bins=(nbins2(y), nbins2(x)),normed=True)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    levels = (0.5,0.4,0.3,0.2,0.15,0.1,0.075,0.06,0.05,0.035,0.02)
    cset = plt.contour(H, levels, origin="lower",linewidths=(1.9, 1.6, 1.5, 1.4),colors="black",extent=extent)
    plt.clabel(cset, inline=1, fontsize=10, fmt="%0.2f")
    for c in cset.collections:
        c.set_linestyle("solid")

class FitParameters:

    exclude=['SIMPLE',
             'BITPIX',
             'NAXIS','NAXIS1','NAXIS2',
             'EXTEND',
             'COMMENT',
             'DIM',
             'H_MIN_X','H_MAX_X','BINSIZE_X',
             'H_MIN_Y','H_MAX_Y','BINSIZE_Y',
             'FITSTAT']

    def __init__(self,phdr,fitted):
        self.initial=dict()
        for key in phdr.keys():
            if not (key in FitParameters.exclude):
                self.initial[key] = phdr[key]

    def info(self):
        print "Best Fit Parameters"
        keyPrint(self.initial)

class FitImage:
    
    colorMaps={'Residual':cm.Reds,
              'model':cm.Greys,
              'observation':cm.Greens}

    def __init__(self,hdus,extension):
        self.img=hdus[extension].data
        self.phdu=hdus[0].header

        self.xbinsize=self.phdu['BINSIZE_X']
        self.ybinsize=self.phdu['BINSIZE_Y']
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

    def plot(self):
        plt.imshow(self.img, interpolation='bicubic', cmap=FitImage.colorMaps[self.name], extent=self.extent)
        plt.xlabel("AXIS 1")
        plt.ylabel("AXIS 2")
        plt.title(self.name)

    def show(self):
        self.plot()
        plt.show()

    def save(self,imgfile):
        plt.figure()
        self.plot()
        savefig(imgfile)

class NumberCounts:

    def __init__(self,hdus):
        self.name=""

class MCMCInfo:
    def __init__(self,hdus):
        self.Parameters=dict()
        self.Rvalues=dict()
        self.ChiSquares=list()
        self.Acceptance=list()
        self.ChainNum=0

        params=hdus[4].data
        pnames=list()
        for col in params.names:
            if "0" in col:
                if not (("CHISQ" in col) | ("ACPT" in col)):
                    ptemp=col.split('0')[0]
                    pnames.append(ptemp)
            if "CHISQ" in col:
                self.ChainNum+=1
                self.ChiSquares.append(params[col])
            elif "ACPT" in col:
                self.Acceptance.append(params[col])                

        chdu=hdus[5]
        rnames=chdu.data.names
        for i in range(0,len(rnames)):
            rs=chdu.data[rnames[i]]
            gpts=numpy.where(rs > 0)
            self.Rvalues[pnames[i]]=rs[gpts]

        pnames=numpy.unique(pnames)
        for i in range(0,len(pnames)):
            p0=pnames[i]
            if((p0 == "PHI") | (p0 == "L")):
                p0=p0+str(0)
            self.Parameters[p0]=list()
            for i in range(0,self.ChainNum):
                self.Parameters[p0].extend(params[p0+str(i)])
            

    def plotFit(self):
        pnames=self.Parameters.keys()
        dim=len(pnames)
        for i in range(0,len(pnames)):
            for j in range(0,len(pnames)):
                if(i >= j):
                    plt.subplot(dim,dim,i*dim+j+1)
                if(i == j):
                    datapts=self.Parameters[pnames[i]]
                    plt.hist(datapts,
                         bins=nbins(datapts),
                         normed=True,
                         histtype='step',
                         color="black")
                    if(j == 0):
                        plt.ylabel("Relative Probability")
                elif(i > j):
                    contour2d(self.Parameters[pnames[j]],self.Parameters[pnames[i]])
                    if(j == 0):
                        plt.ylabel(pnames[i])
                if(i == (len(pnames)-1)):
                    plt.xlabel(pnames[j])

    def showFit(self):
        self.plotFit()
        plt.show()

    def saveFit(self,filename):
        plt.figure(figsize=(12,12))
        self.plotFit()
        plt.savefig(filename)

    def plotR(self):
        for key in self.Rvalues.keys():
            iterations=numpy.arange(1,len(self.Rvalues[key])+1)
            plt.plot(iterations,self.Rvalues[key],label=key)
        plt.legend(loc='upper right')
        plt.xlabel("Iteration")
        plt.ylabel("R Convergence Criterion")

    def showR(self):
        self.plotR()
        plt.show()
    
    def saveR(self,filename):
        self.plotR()
        plt.savefig(filename)

class OutputFile:

    def __init__(self,filename):
        self.filename=filename
        self.hdus=fits.open(filename)
        self.images=dict()
        if(len(self.hdus) == 3):
            self.fitExt=False
            self.runtype='Simulation Only'
            self.parameters=FitParameters(self.hdus[0].header,False)
        else:
            self.fitExt=True
            self.runtype='Parameter Fitting'
            self.parameters=FitParameters(self.hdus[0].header,True)
            self.MCMC=MCMCInfo(self.hdus)
            for i in range(0,3):
                self.images['temp'] = FitImage(self.hdus,i)
                newkey=self.images['temp'].name
                self.images[newkey] = self.images.pop('temp')

    def type(self):
        return self.runtype

    def fit(self):
        return self.fitExt

    def info(self):
        print self.runtype
        keyPrint(self.images)

    def plotImages(self):
        col=0
        for key in self.images.keys():
            col=col+1
            plt.subplot(1,3,col)
            self.images[key].plot()
            if(col > 1):
                plt.ylabel("")

    def showImages(self):
        self.plotImages()
        plt.show()
