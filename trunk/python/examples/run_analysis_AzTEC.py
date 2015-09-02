#!/usr/bin/env python
# modified for AzTEC run by Jed 

import os
import sys
import datetime
import time
from ModelFile import *

## don't touch these for AzTEC runs on default_templates
templatefile="/Users/Jed/SurveySim/trunk/templates/default_templates_v2.fits" 
obsfile="/Users/Jed/SurveySim/trunk/obs/aztec.fits"



## change these for each new run ##
#####################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

outfile='/Users/Jed/AzTEC/ss_2.0_runs/AzTEC_test/output.fits'
mfile='/Users/Jed/AzTEC/ss_2.0_runs/AzTEC_test/model.fits'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#####################################################################


mod=ModelFile()

mod.axis1="ColorF2F1"
mod.axis2="Flux1"

mod.survey['area']=0.7                 # in deg^2
mod.redshift['min'] = 0.1

# FILTERS #
# -------------------------------------------------------------------------------- # 
mod.filters[0].setID('AzTEC_1.1')        
mod.filters[0].limit=1.3               #0.9 - 1.7 mJy  
mod.filters[0].err=0.4         
mod.filters[0].unit='mJy'
mod.filters[0].serr=0.00


mod.filters[1].setID("Spitzer_MIPS_24um")
mod.filters[1].limit=0.03                   # limit is upper/UDS limit from M+2012 
mod.filters[1].err=-0.011                   # from Ivison et al. 2007 
mod.filters[1].unit='mJy'
mod.filters[1].serr=0.00

                                       # avoid log-error w/ non-zero limits
mod.filters[2].setID("Spitzer_IRAC_4") # according to IRAC instrument handbook
mod.filters[2].limit=0.00001           # http://irsa.ipac.caltech.edu/data/SPITZER/...
mod.filters[2].err=0.0                  # low limit set to detect "non-detections" (very low fluxes)
mod.filters[2].unit='mJy'              # don't really need error --> not even using 8um data
mod.filters[2].serr=0.00
# -------------------------------------------------------------------------------- #
#                                    Parameters                                    #
# -------------------------------------------------------------------------------- #

#mod.setLF("ModifiedSchecter") # straight up didn't work with MS
#mod.setLF('Schecter')         # schecter works, but doesn't change anything 
mod.setLF("DoublePowerLaw")    # I like this the best - used by N+13 in IR-lumfunc 

mod.params['Alpha'].value=0.54 # fixed
mod.params['Alpha'].fixed=0
mod.params['Alpha'].pmin=-3.5
mod.params['Alpha'].pmax=3.5

mod.params['Beta'].value=3.03 # fixed
mod.params['Beta'].fixed=0
mod.params['Beta'].pmin=0.56
mod.params['Beta'].pmax=3.5

mod.params['Phi0'].value=-2.4 # fixed 
mod.params['Phi0'].fixed=1
mod.params['Phi0'].pmin=-2.5
mod.params['Phi0'].pmax=-1.9

mod.params['L0'].value=9.5 # fixed 
mod.params['L0'].fixed=1
mod.params['L0'].pmin=9
mod.params['L0'].pmax=11.0

mod.params['P'].pmax=5 #  
mod.params['P'].pmin=-5
mod.params['P'].value=-3.5
mod.params['P'].fixed=0

mod.params['P2'].pmax=2 # fixed
mod.params['P2'].pmin=-8
mod.params['P2'].value=-2.87
mod.params['P2'].fixed=1

mod.params['Q'].pmax=5 # 
mod.params['Q'].pmin=-5
mod.params['Q'].value=4.5
mod.params['Q'].fixed=0

mod.params['Q2'].pmax=5 # fixed 
mod.params['Q2'].pmin=-5
mod.params['Q2'].value=1.49
mod.params['Q2'].fixed=1

mod.params['zbp'].pmax=4 # fixed 
mod.params['zbp'].pmin=0.5
mod.params['zbq'].pmax=4
mod.params['zbq'].pmin=0.5
mod.params['zbp'].value=2.54
mod.params['zbq'].value=1.64
mod.params['zbp'].fixed=1
mod.params['zbq'].fixed=1

mod.params['zbc'].value=2
mod.params['zbc'].pmin=0
mod.params['zbc'].pmax=4
mod.params['zbc'].fixed=0

mod.params['cexp'].pmin=-3  
mod.params['cexp'].pmax=3 
mod.params['cexp'].value=-2 # from RUNS_3 : AzTEC A
mod.params['cexp'].fixed=0

# doesnt seem like these parameters are affecting redshift
# distribution results 
mod.params['fa0'].value=0.28    # from Johnson et al. 2013
mod.params['fa0'].pmin=0.14        # mimimum from J+2013
mod.params['fa0'].pmax=0.5   # maximum from M+2012 fraction of mm-sources w/ 1.4 GHz data 
mod.params['fa0'].fixed=1

mod.params['t1'].value= 1.9   # from previous runs 
mod.params['t1'].pmin = -2
mod.params['t1'].pmax =  2
mod.params['t1'].fixed= 0

mod.params['t2'].value= 0.7392 # from previous runs 
mod.params['t2'].pmin = -2
mod.params['t2'].pmax =  2
mod.params['t2'].fixed= 0

mod.params['zbt'].value= 1.33 # from previous runs 
mod.params['zbt'].pmin = 0
mod.params['zbt'].pmax = 4
mod.params['zbt'].fixed= 0


# -------------------------------------------------------------------------------- #
# MCMC params #  -  not touching any of these
# -------------------------------------------------------------------------------- # 


mod.settings['verbosity']=3
mod.annealing['temp']=.02
mod.annealing['learningRate']=0.1

mod.convergence['CI']=0.90

# -------------------------------------------------------------------------------- #



# run surveysim 
mod.write(mfile)
mod.info()

os.system('SurveySim '+mfile+' '+templatefile+' '+obsfile+' -o '+outfile)

