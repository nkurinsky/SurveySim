#!/usr/bin/env python
# modified for AzTEC run by Jed 

import os
import sys
import datetime
import time
from ModelFile import *

mfile="../../Examples/AzTEC/LHE_B_model.fits"
templatefile="../templates/default_templates.fits" # make sure to follow path
obsfile="../obs/LHE_complete_c.fits"
outfile="../../Examples/AzTEC/LHE_B_output.fits"

## Info on this specific run 
## ------------------------
## run #: 4
## specs: LHE, sources with complete 1.1 mm and 24 um data  
## notes: parameters are fitted
##        no color evolution, alpha and beta are fixed



mod=ModelFile()

mod.axis1="ColorF1F2"
mod.axis2="Flux1"

mod.survey['area']=0.7 # in deg^2

mod.filters[0].setID('AzTEC_1.1')        
mod.filters[0].limit=1.73               # mJy lower flux limit 
mod.filters[0].err=1.2165             # characteristic flux measurement error

mod.filters[1].setID("Spitzer_MIPS_24um")
mod.filters[1].limit=0.021
mod.filters[1].err=0.012611     # as calculated in LHE_analysis 

# avoid log-error w/ non-zero limits
mod.filters[2].setID("Spitzer_IRAC_4") # according to IRAC instrument handbook
mod.filters[2].limit=0.001                # http://irsa.ipac.caltech.edu/data/SPITZER/...
mod.filters[2].err=0.001

#mod.setLF("ModifiedSchecter")
mod.setLF("DoublePowerLaw")

mod.params['Alpha'].value=1.15 # fixed
mod.params['Alpha'].fixed=1

mod.params['Beta'].value=0.5 # fixed
mod.params['Beta'].fixed=1

mod.params['Phi0'].value=-2.29 # fitted 
mod.params['Phi0'].fixed=0
mod.params['Phi0'].pmin=-2.5
mod.params['Phi0'].pmax=-1.9

mod.params['L0'].value=10.12 # fitted 
mod.params['L0'].fixed=0
mod.params['L0'].pmin=10.0
mod.params['L0'].pmax=11.0

mod.params['P'].pmax=2 # fitted 
mod.params['P'].pmin=-10
mod.params['P'].value=-0.57
mod.params['P'].fixed=0

mod.params['P2'].pmax=2 # fitted
mod.params['P2'].pmin=-10
mod.params['P2'].value=-3.92
mod.params['P2'].fixed=0

mod.params['Q'].pmax=5 # fitted 
mod.params['Q'].pmin=-5
mod.params['Q'].value=3.55
mod.params['Q'].fixed=0 

mod.params['Q2'].pmax=5 # fitted 
mod.params['Q2'].pmin=-5
mod.params['Q2'].value=1.62
mod.params['Q2'].fixed=0

mod.params['zbp'].pmax=4 # fitted 
mod.params['zbp'].pmin=0.5
mod.params['zbq'].pmax=4
mod.params['zbq'].pmin=0.5
mod.params['zbp'].value=1.1
mod.params['zbq'].value=1.85
mod.params['zbp'].fixed=0
mod.params['zbq'].fixed=0

mod.params['cexp'].pmin=0.0 # fixed - no color evolution 
mod.params['cexp'].pmax=4.0 
mod.params['cexp'].value=0
mod.params['cexp'].fixed=1

mod.settings['verbosity']=3
mod.annealing['temp']=.02
mod.annealing['learningRate']=0.1

mod.convergence['CI']=0.90

mod.write(mfile)
mod.info()

os.system('SurveySim '+mfile+' '+templatefile+' '+obsfile+' -o '+outfile)
