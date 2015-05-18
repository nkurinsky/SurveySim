#!/usr/bin/env python

import os
import sys
import datetime
import time
from ModelFile import *

mfile="spire_model.fits"
templatefile="../trunk/templates/Kirkpatrick2015_templates.fits"
obsfile="../trunk/obs/spire.fits"
outfile="spire_output.fits"

mod=ModelFile()

mod.axis1="ColorF2F3"

mod.filters[0].setID("SPIRE_250")
mod.filters[1].setID("SPIRE_350")
mod.filters[2].setID("SPIRE_500")

mod.params['Alpha'].value=3
mod.params['Alpha'].fixed=0
mod.params['Beta'].fixed=0
mod.params['Phi0'].fixed=0
mod.params['Phi0'].pmin=-2.5
mod.params['Phi0'].pmax=-1.9
mod.params['L0'].fixed=0
mod.params['L0'].pmin=9
mod.params['L0'].pmax=11

mod.params['P'].pmax=2
mod.params['P'].pmin=-7
mod.params['P2'].pmax=1
mod.params['P2'].pmin=-7
mod.params['P2'].value=0
mod.params['P2'].fixed=1

mod.params['Q'].pmax=8
mod.params['Q'].pmin=0
mod.params['Q2'].pmax=5
mod.params['Q2'].pmin=-2
mod.params['Q2'].value=0
mod.params['Q2'].fixed=1

mod.params['zbp'].pmax=4
mod.params['zbp'].pmin=0
mod.params['zbq'].pmax=4
mod.params['zbq'].pmin=0
mod.params['zbp'].value=2
mod.params['zbq'].value=2
mod.params['zbp'].fixed=0
mod.params['zbq'].fixed=0

mod.params['cexp'].pmin=0
mod.params['cexp'].pmax=4
mod.params['cexp'].value=2
mod.params['cexp'].fixed=0

mod.settings['verbosity']=3
mod.annealing['temp']=.05
mod.annealing['learningRate']=2

mod.convergence['CI']=0.90

mod.write(mfile)
mod.info()

os.system('SurveySim '+mfile+' '+templatefile+' '+obsfile+' -o '+outfile)
