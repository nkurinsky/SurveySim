import pyfits as fits
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
sns.set(style='white')#,font='serif')#,font_scale=1.5,palette='Set2')
#sns.despine()

plt.rc('text', usetex=True)

params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{sfmath}']}
plt.rcParams.update(params)

from ModelFile import keyPrint

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
        else:
            self.name=hdus[extension].header['EXTNAME']


    def __str__(self):
        return "Image: "+self.name

    def plot(self,interpolation='nearest',cmap=cm.Greys,labelCbar=True):
        plt.imshow(self.img, interpolation=interpolation, cmap=cmap, extent=self.extent, aspect='auto')
        cbar=plt.colorbar()
        if(labelCbar):
            cbar.set_label("Normalized Density")
        plt.contour(self.img, interpolation=interpolation, colors='grey', extent=self.extent, origin='image', aspect='auto')
        plt.xlabel("AXIS 1")
        plt.ylabel("AXIS 2")
        plt.title(self.name)

    def show(self):
        self.plot()
        plt.show(block=True)

    def save(self,imgfile):
        plt.figure()
        self.plot()
        savefig(imgfile)

class NumberCounts:

    def __init__(self,hdus,fitted):
        self.chains=dict()
        self.fitted=fitted
        if(fitted):
            chainExt=7
            dists=hdus[3].data
            self.best=dict()
            self.bins=[dists['s1'],dists['s2'],dists['s3']]
            self.best['observed']=[dists['obs_dnds1'],dists['obs_dnds2'],dists['obs_dnds3']]
            self.best['modeled']=[dists['mod_dnds1'],dists['mod_dnds2'],dists['mod_dnds3']]
            self.chains['chisq']=hdus[chainExt].data['chisq']
        else:
            chainExt=2
            dists=hdus[1].data
            self.single=dict()
            self.bins=[dists['s1'],dists['s2'],dists['s3']]
            self.single['modeled']=[dists['mod_dnds1'],dists['mod_dnds2'],dists['mod_dnds3']]

        self.chains['counts']=[hdus[chainExt].data['dnds1'],
                               hdus[chainExt].data['dnds2'],
                               hdus[chainExt].data['dnds3']]

    def info(self):
        print level1,"Number Counts"
        print level2,"Bins (bins)"
        if(self.fitted):
            print level2,"Fitted Counts (best)"
            print level3,self.best.keys()
        else:
            print level2,"Single Run Counts (single)"
            print level3,self.single.keys()
        print level2,"Simulated Count Chains (chains)"
        print level3,self.chains.keys()

    def plot(self,band,showLegend=True):
        bi=band-1
        if(self.fitted):
            chi=self.chains['chisq']
            chimed=np.median(chi)
            gpts=np.where(chi < chimed)
            chainCounts=self.chains['counts'][bi][gpts]
        else:
            chainCounts=self.chains['counts'][bi]
            
        bins=self.bins[bi]
        medians=np.zeros(len(bins))
        errlow=np.zeros(len(bins))
        errhigh=np.zeros(len(bins))
        for i in range(0,len(chainCounts[0])):
            link=chainCounts[:,i]
            medians[i]=np.median(link)
            sortpts=np.argsort(link)
            lowind=int(0.159*len(sortpts))
            highind=int((1.0-0.159)*len(sortpts))
            errlow[i]=medians[i]-link[sortpts[lowind]]
            errhigh[i]=link[sortpts[highind]]-medians[i]

        if(self.fitted):            
            model=self.best['modeled'][bi]
            obs=self.best['observed'][bi]
            opts=np.where(obs > 1e-10)
            pts=np.where(model > 1e-10)
            plt.plot(bins[opts],obs[opts],'o',label="Observation",color='black')
            plt.plot(bins[pts],model[pts],'d',label="Best",color='black')
        else:
            model=self.single['modeled'][bi]
            pts=np.where(model > 1e-10)
            plt.plot(bins[pts],model[pts],'d',label="Best",color='black')
        
        mpts=np.where(medians > 1e-10)
        #sns.set(font='serif',font_scale=1.5,palette='Set2')
        sns.set(style='white',font_scale=1.5,palette='Set2')
        #sns.despine()

        plt.plot(bins[mpts],medians[mpts],'s',label="Median",color='black')
        plt.errorbar(bins[mpts],medians[mpts],yerr=[errlow[mpts],errhigh[mpts]],linestyle="None",color='gray')
        #this needs to be made non-run specific perhaps by repeating the band keywords from the model into this 
        plt.title("Counts Band "+str(band))

        plt.xlabel(r'Flux [mJy]')
        plt.ylabel(r'$S^{2.5} dN/dSd\Omega*[Jy^{1.5} ster^{-1}]$')
        plt.xscale('log')
        plt.yscale('log')

        if(showLegend):
            plt.legend(loc='lower right')

    def show(self,band):
        self.plot(band)
        plt.show(block=True)

    def save(self,band,filename):
        self.plot(band)
        plt.savefig(filename)

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

        chdu=hdus[5]
        rnames=chdu.data.names
        for i in range(0,len(rnames)):
            rs=chdu.data[rnames[i]]
            gpts=np.where(rs > 0)
            self.Rvalues[pnames[i]]=rs[gpts]

        pnames=np.unique(pnames)
        for i in range(0,len(pnames)):
            p0=pnames[i]
            if((p0 == "PHI") | (p0 == "L")):
                p0=p0+str(0)
            self.Parameters[p0]=list()
            for i in range(0,self.ChainNum):
                gpts=np.where(self.ChiSquares[i] < np.median(self.ChiSquares))
                self.Parameters[p0].extend(params[p0+str(i)][gpts])

        for k,v in self.Parameters.iteritems():
            self.BestFit[k]={'mean':np.mean(self.Parameters[k]),
                             'median':np.median(self.Parameters[k]),
                             'sigma':std(self.Parameters[k])}

    def info(self):
        print level2,"Fitted Parameters (BestFit)"
        for k,v in self.BestFit.iteritems(): print level3,k,'-\t',v
        print level2,"Number of Chains:",self.ChainNum
        print level2,"Available Variable Chains:"
        print level3,"Parameters"
        print level3,"Acceptance"
        print level3,"ChiSquares"
        print level3,"Convergence"

    def plotFit(self):
        pnames=self.Parameters.keys()    
        #want to re-order these as they do not appear in the most sensible order
        pnames_default=['PHI0','L0','BETA','ALPHA','P','Q','ZBP','ZBQ','P2','Q2','CEXP','ZBC']
        pnames_sorted=[]
        for p in pnames_default:
            if p in pnames:
                pnames_sorted.append(p)
        pnames=pnames_sorted
        dim=len(pnames)
        for i in range(0,dim):
            for j in range(0,dim):
                if(i >= j):
                    ax=plt.subplot(dim,dim,i*dim+j+1)
                    plt.subplots_adjust(hspace = .05, wspace=0.05)
                if(i == j):
                    datapts=self.Parameters[pnames[i]]
                    #plt.hist(datapts,
                    #     bins=nbins(datapts),
                    #     normed=True,
                    #     histtype='step',
                    #     color="black")
                    datapts=np.array(datapts)
                    datapts.astype(float)
                    sns.despine()
                    sns.kdeplot(datapts)
                    if(j == 0):
                        plt.ylabel("Probability")
                    else:
                        plt.setp( plt.gca().get_yticklabels(), visible=False)
                elif(i > j):
                    x=self.Parameters[pnames[j]]
                    y=self.Parameters[pnames[i]]
                    x=np.array(x)
                    y=np.array(y)
                    x.astype(float)
                    y.astype(float)
                    #MCMCcontour(self.Parameters[pnames[j]],self.Parameters[pnames[i]])
                    sns.kdeplot(x,y,n_levels=10,cmap="Purples_d")
                    start, end = ax.get_xlim()
                    stepsize=(end-start)/5. #force to only show a maximum of 5 tickmarks
                    if(stepsize >= 1):
                        stepsize=int(stepsize+0.5)
                    if((stepsize > 0.3) & (stepsize < 1)):
                        stepsize=0.5
                    xticks=np.arange(start, end, stepsize)
                    xticks=np.around(xticks)
                    ax.xaxis.set_ticks(xticks)
                    
                    start, end = ax.get_ylim()
                    stepsize=(end-start)/5. #force to only show a maximum of 5 tickmarks
                    yticks=np.around(np.arange(start, end, stepsize)+0.5,1)
                    ax.yaxis.set_ticks(yticks)
                    if(j == 0):
                        tmp=pnames[i]
                        if(tmp == 'ALPHA'):
                            plt.ylabel(r'$\alpha$')
                        if(tmp == 'BETA'):
                            plt.ylabel(r'$\beta$')
                        if(tmp == 'PHI0'):
                             plt.ylabel(r'$\Phi_{*,0}$')
                        if(tmp == 'L0'):
                             plt.ylabel(r'$L_{*,0}$')
                        if(tmp == 'P'):
                            plt.ylabel(r'$p_1$')
                        if(tmp == 'Q'):
                            plt.ylabel(r'$q_1$')
                        if(tmp == 'P2'):
                            plt.ylabel(r'$p_2$')
                        if(tmp == 'Q2'):
                            plt.ylabel(r'$q_2$')                           
                        if(tmp == 'ZBP'):
                            plt.ylabel(r'$z_{break,p}$')
                        if(tmp == 'ZBQ'):
                            plt.ylabel(r'$z_{break,q}$')
                        if(tmp == 'ZBC'):
                            plt.ylabel(r'$z_{break,c}$')
                        if(tmp == 'CEXP'):
                            plt.ylabel(r'$c_{exp}$')
                    else:
                        plt.setp( plt.gca().get_yticklabels(), visible=False)
                if(i == (len(pnames)-1)):
                    tmp=pnames[j]
                    if(tmp == 'ALPHA'):
                        plt.xlabel(r'$\alpha$')
                    if(tmp == 'BETA'):
                        plt.xlabel(r'$\beta$')
                    if(tmp == 'PHI0'):
                        plt.xlabel(r'$\Phi_{*,0}$')
                    if(tmp == 'L0'):
                        plt.xlabel(r'$L_{*,0}$')
                    if(tmp == 'P'):
                        plt.xlabel(r'$p_1$')
                    if(tmp == 'Q'):
                        plt.xlabel(r'$q_1$')
                    if(tmp == 'P2'):
                        plt.xlabel(r'$p_2$')
                    if(tmp == 'Q2'):
                        plt.xlabel(r'$q_2$')
                    if(tmp == 'ZBP'):
                        plt.xlabel(r'$z_{break,p}$')
                    if(tmp == 'ZBQ'):
                        plt.xlabel(r'$z_{break,q}$')
                    if(tmp == 'ZBC'):
                        plt.xlabel(r'$z_{break,c}$')
                    if(tmp == 'CEXP'):
                        plt.xlabel(r'$c_{exp}$')
                else:
                    plt.setp( plt.gca().get_xticklabels(), visible=False)
        #add legend with best-fit parameter values
        ax.annotate('Best-fit values:', xy=(.7, .8),  xycoords='figure fraction',
                horizontalalignment='center', verticalalignment='center',size=20)
        yoffset=0.03
        for k,v in self.BestFit.iteritems(): 
        #    print level3,k,'-\t',v
             #ax.annotate(k,'=',xy=(.7,.8-yoffset),xycoords='figure fraction',
             #   horizontalalignment='center', verticalalignment='center',size=20)
             yoffset=yoffset+0.03
        #print self.BestFit(pnames[0])


    def showFit(self,block=True):
        my_dpi=116
        maxsize=maxWindowSize()
        figsize=(maxsize[0]/float(my_dpi)-0.25,maxsize[1]/float(my_dpi)-0.5)
        #plt.figure(figsize=figsize)
        fig=plt.figure(figsize=[12,12])
        self.plotFit()
        plt.show(block=block)

    def saveFit(self,filename):
        plt.figure(figsize=(12,12))
        self.plotFit()
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
        self.filename=filename
        self.hdus=fits.open(filename)
        self.images=dict()
        if(len(self.hdus) == 3):
            self.fitExt=False
            self.runtype='Simulation Only'
            self.simInfo=FitInfo(self.hdus,False)
            self.counts=NumberCounts(self.hdus,False)
        else:
            self.fitExt=True
            self.runtype='Parameter Fitting'
            self.simInfo=FitInfo(self.hdus,True)
            self.counts=NumberCounts(self.hdus,True)
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
        self.counts.info()

    def plotImages(self,xrange=None,yrange=None):
        col=0
        for key in self.images.keys():
            col=col+1
            plt.subplot(1,3,col)
            self.images[key].plot(labelCbar=False)
            if(col > 1):
                plt.ylabel("")
            if(xrange != None):
                plt.xlim(xrange[0],xrange[1])
            if(yrange != None):
                plt.ylim(yrange[0],yrange[1])
        #plt.gcf().set_tight_layout(True)
        plt.tight_layout(True)

    def showImages(self,block=True,xrange=None,yrange=None):
        plt.figure(figsize=(12,4))
        self.plotImages(xrange=xrange,yrange=yrange)
        plt.show(block=block)

    def saveImages(self,filename,xrange=None,yrange=None):
        plt.figure(figsize=(12,4))
        self.plotImages(xrange=xrange,yrange=yrange)
        plt.savefig(filename)

    def plotCounts(self):
        for i in range(1,4):
            plt.subplot(1,3,i)
            self.counts.plot(i)
            if(i > 1):
                plt.ylabel("")

    def showCounts(self,block=True):
        plt.figure(figsize=(16,6))
        self.plotCounts()
        plt.show(block=block)

    def saveCounts(self,filename):
        plt.figure(figsize=(12,5))
        self.plotCounts()
        plt.savefig(filename)

    def plot(self):
        if(self.fit()):
            rows=3
        else:
            rows=2
        for i in range(1,4):
            plt.subplot(rows,3,i)
            self.counts.plot(i,showLegend=False)
            if(i > 1):
                plt.ylabel("")
        
        plt.subplot(rows,3,4)
        self.simInfo.plotRedshift()
        plt.subplot(rows,3,5)
        self.simInfo.plotLuminosity()
        plt.subplot(rows,3,6)
        self.simInfo.plotSEDType()

        col=0
        for key in self.images.keys():
            col=col+1
            plt.subplot(rows,3,col+6)
            self.images[key].plot(labelCbar=False)
            if(col > 1):
                plt.ylabel("")
        #plt.gcf().set_tight_layout(True)
        plt.tight_layout(True)

    def show(self,block=True):
        my_dpi=116 
        maxsize=maxWindowSize()
        figsize=(maxsize[0]/float(my_dpi)-0.25,maxsize[1]/float(my_dpi)-0.5)
        plt.figure(figsize=figsize, dpi=my_dpi)
        self.plot()
        plt.show(block=block)
