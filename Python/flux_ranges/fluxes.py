import pyfits
import numpy
from pylab import *

from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u

hdus = pyfits.open("/usr/local/surveysim/templates/sf_templates.fits")
ftable = hdus[0].data

wave = ftable[0]*u.micrometer
lum = ftable[14]*u.W/u.Hz

color_evolution_options = numpy.arange(0,3.01,0.25)
redshifts = numpy.arange(0.0,4.01,0.1)
bands = [250.0,350.0,500.0]

for bandi in bands:
    for modi in range(1,15):
        clf()
        lum = ftable[modi]*u.W/u.Hz
        for i in redshifts:
            emitted = bandi*u.micrometer/(1.0+i)
            gpts = where(wave > emitted)
            mylum = lum[gpts[0][0]]
            print(wave[gpts[0][0]])
            dist = cosmo.luminosity_distance(i)
            dist = dist.to(u.m)
            print(mylum)
            flux = (1+i)*mylum/(4*3.1415*pow(dist,2))
            flux = flux.to(u.mJy)
            print(flux)
            for j in color_evolution_options:
                if(i < 2):
                    scatter(i,log10(flux.value*pow((1+i),j)),j*2)
                else:
                    scatter(i,log10(flux.value*pow(3,j)),j*2)
                    
        ylabel("Log (mJy)")
        xlabel("Redshift")
        xlim(0,4.1)
        ylim(1.0,4.0)
        savefig("fluxes_model"+str(modi)+"_band"+str(bandi)+".eps")
