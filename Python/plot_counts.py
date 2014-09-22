import pyfits as pf
import numpy as np
import matplotlib.pyplot as p

def plot_counts(filename,bandnum,alpha=0.159):

    hdus = pf.open(filename)
    dists = hdus[3].data
    count_dists = hdus[6].data

    plusfrac = (1.0-alpha)
    minusfrac = alpha

    b = str(bandnum)
    if((bandnum > 0) and (bandnum < 4)):
        x = dists['s'+b]
        ybest = dists['mod_dnds'+b]
        yobs = dists['obs_dnds'+b]
    else: 
        print "Invalid band number: "+b
        return(1)

    c = count_dists["dnds"+b]
    chis = count_dists["chisq"]
    csize = len(c[0])
    cmean = np.zeros(csize)
    cplus = np.zeros(csize)
    cminus = np.zeros(csize)

    goodpts = np.where(chis < np.median(chis))[0]

    for i in range(0,csize):
        dnds = c[goodpts,i]
        dnds = np.sort(dnds)
        gpts = np.where(dnds > 0)
        dnds = dnds[gpts[0]]
        if(len(dnds) > 0):
            pi = int(plusfrac*len(dnds))
            mi = int(minusfrac*len(dnds))
            cmean[i] = np.median(dnds)
            cplus[i] = dnds[pi]
            cminus[i] = dnds[mi]
    
    bestpts = np.where(ybest > 0)
    xbest = x[bestpts[0]]/1e3
    ybest = ybest[bestpts[0]]
    
    obspts = np.where(yobs > 0)
    xobs = x[obspts[0]]/1e3
    yobs = yobs[obspts[0]]
    
    minpts = np.where(cminus > 0)
    xstats = x[minpts[0]]/1e3
    yminus = cminus[minpts[0]]
    yplus = cplus[minpts[0]]
    ymean = cmean[minpts[0]]
    
    p.clf()
    p.xscale('log')
    p.yscale('log')
    #p.loglog(xstats,ymean,ls='-',label="median",color="black")
    #p.loglog(xstats,yplus,ls='--',label="$1\sigma$",color="black")
    #p.loglog(xstats,yminus,ls='--',color="black")
    p.scatter(xobs,yobs,label="observed",color="blue")
    p.scatter(xbest,ybest,label="model",color="green")
   
    p.legend()
    
    p.xlabel("F"+b+" [Jy]")
    p.ylabel("$(dN/dS)S^{2.5} [gal\  ster^{-1} J^{1.5}]$")
    return(p)

def all():
    for i in range(1,4):
        p = plot_counts("output.fits",i)
        p.savefig("band"+str(i)+".png")
