#!/usr/bin/env python

import os
import sys
import datetime
import time
from ModelFile import *

mfile="spire_model.fits"
templatefile="../templates/Kirkpatrick2015_templates.fits"
obsfile="../obs/spire_cdfs-swire.fits"
outfile="spire_output.fits"

mod=ModelFile()

mod.axis1="ColorF2F3"

mod.survey['area']=11.1

mod.filters[0].setID("SPIRE_250")
mod.filters[0].limit=2
mod.filters[0].err=2.5

mod.filters[1].setID("SPIRE_350")
mod.filters[1].limit=4
mod.filters[1].err=2.1

mod.filters[2].setID("SPIRE_500")
mod.filters[2].limit=0.1
mod.filters[2].err=3

#mod.setLF("ModifiedSchecter")
mod.setLF("DoublePowerLaw")

mod.params['Alpha'].value=3.03
mod.params['Alpha'].fixed=1

mod.params['Beta'].value=0.52
mod.params['Beta'].fixed=1

mod.params['Phi0'].value=-2.29
mod.params['Phi0'].fixed=0
mod.params['Phi0'].pmin=-2.5
mod.params['Phi0'].pmax=-1.9

mod.params['L0'].value=10.6
mod.params['L0'].fixed=0
mod.params['L0'].pmin=10.0
mod.params['L0'].pmax=11.0

mod.params['P'].pmax=-5
mod.params['P'].pmin=-10
mod.params['P2'].pmax=2
mod.params['P2'].pmin=-10
mod.params['P2'].value=0
mod.params['P2'].fixed=0

mod.params['Q'].pmax=5
mod.params['Q'].pmin=-5
mod.params['Q2'].pmax=5
mod.params['Q2'].pmin=-5
mod.params['Q2'].value=0
mod.params['Q2'].fixed=0

mod.params['zbp'].pmax=4
mod.params['zbp'].pmin=0.5
mod.params['zbq'].pmax=4
mod.params['zbq'].pmin=0.5
mod.params['zbp'].value=1.1
mod.params['zbq'].value=1.85
mod.params['zbp'].fixed=0
mod.params['zbq'].fixed=0

mod.params['cexp'].pmin=0.0
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
