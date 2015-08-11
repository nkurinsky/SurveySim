#!/usr/bin/env python

import os
import sys
import datetime
import time
from ModelFile import *

ssdir="/usr/local/surveysim/"

mfile="spire_model.fits"
templatefile=ssdir+"templates/Kirkpatrick2015_templates.fits"
obsfile=ssdir+"obs/spire_cdfs-swire.fits"
outfile="spire_output.fits"

mod=ModelFile()

mod.axis1="ColorF1F2"
mod.axis2="Flux1"

mod.survey['area']=11.1

mod.filters[0].setID("SPIRE_250")
mod.filters[0].limit=12.7
mod.filters[0].err=2.5

mod.filters[1].setID("SPIRE_350")
mod.filters[1].limit=0.01
mod.filters[1].err=2.1
#mod.filters[1].serr=1.0

mod.filters[2].setID("SPIRE_500")
mod.filters[2].limit=0.01
mod.filters[2].err=3.0

#mod.setLF("ModifiedSchecter")
mod.setLF("DoublePowerLaw")

mod.params['Alpha'].value=0.5
mod.params['Alpha'].pmin=0
mod.params['Alpha'].pmax=2.0
mod.params['Alpha'].fixed=0

mod.params['Beta'].value=3.0
mod.params['Beta'].pmin=2.0
mod.params['Beta'].pmax=6.0
mod.params['Beta'].fixed=0

mod.params['Phi0'].value=-3.1
mod.params['Phi0'].fixed=0
mod.params['Phi0'].pmin=-5.0
mod.params['Phi0'].pmax=-2.0

mod.params['L0'].value=11.2
mod.params['L0'].fixed=0
mod.params['L0'].pmin=10.0
mod.params['L0'].pmax=12.0

mod.params['P'].value=-0.57
mod.params['P'].pmax=-1
mod.params['P'].pmin=-8
mod.params['P'].fixed=0
mod.params['P2'].pmax=2
mod.params['P2'].pmin=-10
mod.params['P2'].value=-3.92
mod.params['P2'].fixed=0

mod.params['Q'].value=3.55
mod.params['Q'].pmax=8
mod.params['Q'].pmin=-2
mod.params['Q'].fixed=0
mod.params['Q2'].pmax=5
mod.params['Q2'].pmin=-5
mod.params['Q2'].value=1.62
mod.params['Q2'].fixed=0

mod.params['zbp'].pmax=4
mod.params['zbp'].pmin=0
mod.params['zbq'].pmax=4
mod.params['zbq'].pmin=0
mod.params['zbp'].value=1.1
mod.params['zbq'].value=1.85
mod.params['zbp'].fixed=0
mod.params['zbq'].fixed=0

mod.settings['verbosity']=3
mod.annealing['temp']=.03
mod.annealing['learningRate']=0.4

mod.convergence['CI']=0.95

simname="spire_DPL_Br"
mfile=simname+"_model.fits"
outfile=simname+"_output.fits"

mod.params['cexp'].value=0
mod.params['cexp'].fixed=1
mod.params['zbc'].value=0
mod.params['zbc'].fixed=1
mod.write(mfile)
mod.info()

os.system('SurveySim '+mfile+' '+templatefile+' '+obsfile+' -o '+outfile)
