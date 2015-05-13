#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("/usr/local/surveysim/python")
from ModelFile import *

mod=ModelFile()
mod.filters[0].setID("SPIRE_250")
mod.filters[1].setID("SPIRE_350")
mod.filters[2].setID("SPIRE_500")

mod.params['P'].pmax=2
mod.params['P'].pmin=-6
mod.params['P2'].pmax=2
mod.params['P2'].pmin=-6

mod.params['Q'].pmax=8
mod.params['Q'].pmin=0
mod.params['Q2'].pmax=6
mod.params['Q2'].pmin=-2

mod.params['zbp'].fixed=0
mod.params['zbq'].fixed=0

mod.settings['verbosity']=3
mod.annealing['tmax']=20
mod.annealing['tscale']=0.001

mod.convergence['r_max']=1.10
mod.convergence['CI']=0.1

mod.write('spire_model.fits')
mod.info()

templatefile="../trunk/templates/Kirkpatrick2015_templates.fits"
modfile="spire_model.fits"
obsfile="../trunk/obs/spire.fits"
outfile="spire_output.fits"

os.system('SurveySim '+modfile+' '+templatefile+' '+obsfile+' -o '+outfile)
