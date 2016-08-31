from astropy.io import fits
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import numpy as np
from numpy import std
from math import ceil
import os
import Tkinter as tk

#Prelims to make prettier colors
import seaborn as sns
sns.set(style='white')

plt.rc('text', usetex=True)

params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{sfmath}']}
plt.rcParams.update(params)

from SurveySim.ModelFile import keyPrint

level1="  "
level2="    "
level3="      "
level4="        "

def binSize (dat): 
    return 3.49*std(dat)/pow(len(dat),0.333)

def nbins (dat): 
    return int(ceil((max(dat)-min(dat))/binSize(dat)))

def nbins2d (dat): 
    return int(ceil((max(dat)-min(dat))/binSize(dat)/2))

def MCMCcontour(x,y):
    H, xedges, yedges = np.histogram2d(y, x, bins=(nbins2d(y), nbins2d(x)),normed=True)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    hmax=np.max(H)
    levels=[0.95*hmax,0.9*hmax,0.8*hmax,0.7*hmax,0.6*hmax,0.5*hmax]
    cset = plt.contour(H, origin="lower",colors="black",extent=extent,levels=levels)

def maxWindowSize():
    root=tk.Tk()
    root.state('withdrawn')
    size=root.maxsize()
    return size

def axisLabel(pname):
    if(pname == 'ALPHA'):
        return r'$\alpha$'
    if(pname == 'BETA'):
        return r'$\beta$'
    if(pname == 'PHI0'):
        return r'$\Phi_{*,0}$'
    if(pname == 'L0'):
        return r'$L_{*,0}$'
    if(pname == 'P'):
        return r'$p_1$'
    if(pname == 'Q'):
        return r'$q_1$'
    if(pname == 'P2'):
        return r'$p_2$'
    if(pname == 'Q2'):
        return r'$q_2$'                           
    if(pname == 'ZBP'):
        return r'$z_{break,p}$'
    if(pname == 'ZBQ'):
        return r'$z_{break,q}$'
    if(pname == 'ZBC'):
        return r'$z_{break,c}$'
    if(pname == 'CEXP'):
        return r'$c_{exp}$'
    if(pname == 'FA0'):
        return r'$f_{agn,0}$'
    if(pname == 'T1'):
        return r'$t_{agn,1}$'
    if(pname == 'T2'):
        return r'$t_{agn,2}$'
    if(pname == 'ZBT'):
        return r'$z_{break,t}$'
    if(pname == 'FCOMP'):
        return r'$f_{comp}$'
    if(pname == 'FCOLD'):
        return r'$f_{cold}$'
    else:
        return pname
    
class FitInfo:

    exclude=['SIMPLE',
             'BITPIX',
             'NAXIS','NAXIS1','NAXIS2',
             'EXTEND',
             'COMMENT',
             'DIM',
             'H_MIN_X','H_MAX_X','BINSIZE_X',
             'H_MIN_Y','H_MAX_Y','BINSIZE_Y',
             'FITSTAT']

    def __init__(self,hdus,fitted):
        phdr=hdus[0].header
        self.fitstat=phdr['FITSTAT']
        self.initial=dict()
        for key in phdr.keys():
            if not (key in FitInfo.exclude):
                self.initial[key] = phdr[key]
        if(fitted):
            paramExt=3
        else:
            paramExt=1
        dists=hdus[paramExt].data
        self.redshift=dists['z']
        self.lums=dists['Lum']
        if( 'Type' in dists.names): #may come out without type depending on SED Library
            self.sedType=dists['Type']
        else:
            self.sedType=[0,0,0,0]
        self.fluxes=[dists['F1'],dists['F2'],dists['F3']]

    def info(self):
        print level2,"Best Fit ChiSquare: "+str(self.fitstat)
        print level2,"Code Fitted Parameters: (initial)"
        for k,v in self.initial.iteritems(): print level3,k,'-\t',v
        print level2,"Variables:"
        print level3,"Redshift (redshift)"
        print level3,"Luminosity (lums)"
        print level3,"Fluxes (fluxes)",len(self.fluxes),len(self.fluxes[0])

    def plotLuminosity(self):
        plt.hist(self.lums,bins=nbins(self.lums),histtype='step',color="black");
        plt.xlabel('Lum')
        plt.ylabel('N(Lum)')

    def showLuminosity(self):
        self.plotLuminosity()
        plt.show(block=True)
        
    def saveLuminosity(self,filename):
        self.plotLuminosity()
        plt.savefig(filename)

    def plotSEDType(self):
        bins=max(self.sedType)-min(self.sedType)+1
        plt.hist(self.sedType,bins=bins,histtype='step',color="black",log=True);
        plt.xlabel('SED Type')
        plt.ylabel('N(Type)')

    def showSEDType(self):
        self.plotSEDType()
        plt.show(block=True)
        
    def saveSEDType(self,filename):
        self.plotSEDType()
        plt.savefig(filename)

    def plotRedshift(self):
        plt.hist(self.redshift,bins=nbins(self.redshift),histtype='step',color="black");
        plt.xlabel('Redshift')
        plt.ylabel('N(z)')

    def showRedshift(self):
        self.plotRedshift()
        plt.show(block=True)
        
    def saveRedshift(self,filename):
        self.plotRedshift()
        plt.savefig(filename)

    def plot(self):
        plt.subplot(1,3,1)
        self.plotRedshift()
        plt.subplot(1,3,2)
        self.plotLuminosity()
        plt.subplot(1,3,3)
        self.plotSEDType()
        
    def show(self):
        self.plot()
        plt.show(block=True)

    def save(self,filename):
        self.plot()
        plt.savefig(filename)

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
            self.img=hdus[2].data-hdus[1].data
        else:
            self.name=hdus[extension].header['EXTNAME']

        self.axis1='Axis 1'
        self.axis2='Axis 2'

    def __str__(self):
        return "Image: "+self.name

    def plot(self,interpolation='nearest',cmap=None,labelCbar=True,cmax=None):

        tmpimg=np.flipud(self.img)
        if(cmap == None):
            cmap=cm.Greys
            if(self.name == 'Residual'):
                cmap=cm.bwr
            elif(self.name == 'Model Diagnostic'):
                cmap=cm.Blues
            elif(self.name=='Observation Diagnostic'):
                cmap=cm.Reds

        if(cmax != None):
            clim=[0,cmax]
            if(self.name == 'Residual'):
                clim=[-cmax,cmax]           

        plt.imshow(tmpimg, interpolation=interpolation, cmap=cmap, extent=self.extent, aspect='auto')
        if(cmax != None):
            plt.clim(clim)

        cbar=plt.colorbar()
        if(labelCbar):
            cbar.set_label("Normalized Density")
        plt.contour(tmpimg, interpolation=interpolation, colors='grey', extent=self.extent, origin='image', aspect='auto')
        plt.xlabel(self.axis1)
        plt.ylabel(self.axis2)

        if(self.name == 'Residual'):
            gpts=(tmpimg > 0)
            print 'St dev of residual: ',np.std(tmpimg[gpts])
            std_str=np.str(np.std(tmpimg[gpts]))
            plt.title(self.name+': st.dev='+std_str[:4])
        else:
            plt.title(self.name)

    def show(self):
        self.plot()
        plt.show(block=True)

    def save(self,imgfile):
        plt.figure()
        self.plot()
        savefig(imgfile)

class MCMCInfo:
    def __init__(self,hdus):
        self.Parameters=dict()
        self.BestFit=dict()
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

        cmed=np.median(self.ChiSquares[self.ChiSquares < 50])

        chdu=hdus[5]
        rnames=chdu.data.names
        for i in range(0,len(rnames)):
            rs=chdu.data[rnames[i]]
            gpts=np.where(rs > 0)
            self.Rvalues[pnames[i]]=rs[gpts]

        pnames=np.unique(pnames)
        for i in range(0,len(pnames)):
            p0=pnames[i]
            if((p0 == "PHI") | (p0 == "L") | (p0 == "FA")):
                p0=p0+str(0)
            self.Parameters[p0]=list()
            for i in range(0,self.ChainNum):
                gpts=np.where(self.ChiSquares[i] < cmed)
                self.Parameters[p0].extend(params[p0+str(i)][gpts])
        gpts=np.where(np.array(self.Parameters[self.Parameters.keys()[0]]) != 0)
        for k in self.Parameters.keys():
            self.Parameters[k]=(np.array(self.Parameters[k])[gpts[0]]).tolist()

        for k,v in self.Parameters.iteritems():
            self.BestFit[k]={'mean':np.mean(self.Parameters[k]),
                             'median':np.median(self.Parameters[k]),
                             'sigma':std(self.Parameters[k])}

    def info(self):
        print level2,"Fitted Parameters (BestFit)"
        for k,v in self.BestFit.iteritems(): print level3,k,'-\t',v
        print level2,"Number of Chains:",self.ChainNum

    def plotFit(self,mode='all'):
        pnames=self.Parameters.keys()    
        #want to re-order these as they do not appear in the most sensible order

        if(mode == 'all'):
            pnames_default=['PHI0','L0','BETA','ALPHA','P','Q','P2','Q2','ZBP','ZBQ','CEXP','ZBC','FA0','T1','T2','ZBT','FCOMP','FCOLD']
        elif(mode == 'lf'):
            pnames_default=['PHI0','L0','BETA','ALPHA','P','Q','P2','Q2','ZBP','ZBQ','CEXP','ZBC']
        elif(mode == 'sed'):
            pnames_default=['FA0','T1','T2','ZBT','FCOMP','FCOLD']
        elif(mode == 'sed_short'):
            pnames_default=['FA0','FCOMP']

        else:
            raise ValueError("Invalid Plot Mode \""+mode+"\"")

        
        pnames_sorted=[]
        for p in pnames_default:
            if p in pnames:
                pnames_sorted.append(p)
        pnames=pnames_sorted
        if((len(pnames) != len(self.Parameters.keys())) and mode=='all'):
            raise Warning('Not All Parameters Found')

        print 'pnames',pnames
        
        dim=len(pnames)
        for i in range(0,dim):
            for j in range(0,dim):
                if(i >= j):
                    ax=plt.subplot(dim,dim,i*dim+j+1)
                    plt.subplots_adjust(hspace = .05, wspace=0.05)
                if(i == j):
                    datapts=self.Parameters[pnames[i]]
                    start = np.min(datapts)
                    end = np.max(datapts)
                    datapts=np.array(datapts)
                    datapts.astype(float)
                    sns.despine()
                    sns.kdeplot(datapts)
                    plt.xlim(start,end)
                    if(j == 0):
                        plt.ylabel("Probability",fontsize=22)
                    else:
                        plt.setp( plt.gca().get_yticklabels(), fontsize=28, visible=False)
                elif(i > j):
                    x=self.Parameters[pnames[j]]
                    y=self.Parameters[pnames[i]]
                    x=np.array(x)
                    y=np.array(y)
                    x.astype(float)
                    y.astype(float)
                    sns.kdeplot(x,y,n_levels=10,cmap="Purples_d")

                    start = np.min(x)
                    end = np.max(x)
                    stepsize=(end-start)/3. #force to only show a maximum of 4 tickmarks
                    if(stepsize >= 1):
                        stepsize=int(stepsize+0.5)

                    xticks=np.arange(start, end, stepsize)
                    xticks=np.around(xticks,2)
                    ax.xaxis.set_ticks(xticks)
                    plt.xlim(start,end)
                    
                    start=np.min(y)
                    end=np.max(y)
                    stepsize=(end-start)/3. #force to only show a maximum of 4 tickmarks
                    yticks=np.around(np.arange(start, end, stepsize),2)
                    ax.yaxis.set_ticks(yticks)
                    if(j == 0):
                        plt.ylabel(axisLabel(pnames[i]),fontsize=28)
                    else:
                        plt.setp( plt.gca().get_yticklabels(), visible=False)
                    plt.ylim(start,end)
                if(i == (len(pnames)-1)):
                    plt.xlabel(axisLabel(pnames[j]),fontsize=28)
                else:
                    plt.setp( plt.gca().get_xticklabels(), fontsize=28,visible=False)

    def showFit(self,block=True,mode='all'):
        my_dpi=116
        maxsize=maxWindowSize()
        figsize=(maxsize[0]/float(my_dpi)-0.25,maxsize[1]/float(my_dpi)-0.5)
        #plt.figure(figsize=figsize)
        fig=plt.figure(figsize=[12,12])
        self.plotFit(mode=mode)
        plt.show(block=block)

    def saveFit(self,filename,mode='all'):
        plt.figure(figsize=(12,12))
        plt.clf()
        self.plotFit(mode=mode)
        plt.savefig(filename)

    def plotR(self):
        for key in self.Rvalues.keys():
            iterations=np.arange(1,len(self.Rvalues[key])+1)
            plt.plot(iterations,self.Rvalues[key],label=key)
        plt.legend(loc='upper right')
        plt.xlabel("Iteration")
        plt.ylabel("R Convergence Criterion")
        plt.xlim(1,max(iterations))

    def showR(self,block=True):
        plt.figure()
        self.plotR()
        plt.show(block=block)
    
    def saveR(self,filename):
        self.plotR()
        plt.savefig(filename)

    def plotChisq(self):
        colors=['black','cyan','red','green','orange','purple','grey','magenta']
        for i in range(0,len(self.ChiSquares)):
            iterations=np.arange(1,len(self.ChiSquares[i])+1)
            gpts = np.where(self.ChiSquares[i] < 100)
            plt.plot(iterations[gpts],self.ChiSquares[i][gpts],label=str(i),color=colors[i])
            gpts = np.where(self.Acceptance[i] > 0)
            plt.plot(iterations[gpts],self.ChiSquares[i][gpts],'--',color=colors[i])
        plt.legend(loc='upper right',title="Chain Number")
        plt.xlabel("Iteration")
        plt.ylabel("Chi-Square Statistic")
        plt.yscale('log')
        plt.ylim(1e-1,10)
        plt.xlim(1,max(iterations))

    def showChisq(self,block=True):
        plt.figure()
        self.plotChisq()
        plt.show(block=block)
    
    def saveChisq(self,filename):
        self.plotChisq()
        plt.savefig(filename)

    def plotChains(self):
        plt.subplot(2,1,1)
        self.plotR()
        plt.subplot(2,1,2)
        self.plotChisq()

    def showChains(self,block=True):
        plt.figure(figsize=(12,8))
        self.plotChains()
        plt.show(block=block)

class OutputFile:

    def __init__(self,filename):
        fnamecore=filename.split('output.fits')[0]
        self.filename=filename
        self.hdus=fits.open(filename)
        self.images=dict()
        if(len(self.hdus) == 3):
            self.fitExt=False
            self.runtype='Simulation Only'
            self.simInfo=FitInfo(self.hdus,False)
        else:
            self.fitExt=True
            self.runtype='Parameter Fitting'
            self.simInfo=FitInfo(self.hdus,True)
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
        print "File Name:",self.filename
        print level1,self.runtype
        self.simInfo.info()
        print level1,"Diagnostic Histograms"
        keyPrint(self.images)
        if(self.fitExt):
            print level1,"MCMC Info"
            self.MCMC.info()

    def plotImages(self,xrange=None,yrange=None):
        for key in self.images.keys():
            if(key == "Observation Diagnostic"):
                col=1
            if(key == "Residual"):
                col=3
            if(key == "Model Diagnostic"):
                col=2
            plt.subplot(1,3,col)
            self.images[key].plot(labelCbar=False)
            if(col > 1):
                plt.ylabel("")
            if(xrange != None):
                plt.xlim(xrange[0],xrange[1])
            if(yrange != None):
                plt.ylim(yrange[0],yrange[1])

        plt.tight_layout(True)

    def showImages(self,block=True,xrange=None,yrange=None):
        plt.figure(figsize=(12,4))
        self.plotImages(xrange=xrange,yrange=yrange)
        plt.show(block=block)

    def saveImages(self,filename,xrange=None,yrange=None):
        plt.figure(figsize=(12,4))
        self.plotImages(xrange=xrange,yrange=yrange)
        plt.savefig(filename)

    def plot(self):
        if(self.fit()):
            rows=2
        else:
            rows=1
        
        plt.subplot(rows,3,1)
        self.simInfo.plotRedshift()
        plt.subplot(rows,3,2)
        self.simInfo.plotLuminosity()
        plt.subplot(rows,3,3)
        self.simInfo.plotSEDType()

        for key in self.images.keys():
            if(key == "Observation Diagnostic"):
                col=1
            if(key == "Residual"):
                col=3
            if(key == "Model Diagnostic"):
                col=2
            plt.subplot(rows,3,col+3)
            self.images[key].plot(labelCbar=False)
            if(col > 1):
                plt.ylabel("")
        plt.tight_layout(True)

    def show(self,block=True):
        my_dpi=116 
        maxsize=maxWindowSize()
        figsize=(maxsize[0]/float(my_dpi)-0.25,maxsize[1]/float(my_dpi)-0.5)
        plt.figure(figsize=figsize, dpi=my_dpi)
        self.plot()
        plt.show(block=block)
