#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from ModelFile import *

mod=ModelFile()
mod.filters[0].setID("SPIRE_250")
mod.filters[1].setID("SPIRE_350")
mod.filters[2].setID("SPIRE_500")
mod.settings['verbosity']=3

mod.write('spire_model.fits')
mod.info()

templatefile="../trunk/templates/Kirkpatrick2015_templates.fits"
modfile="spire_model.fits"
obsfile="../trunk/obs/spire.fits"
outfile="spire_output.fits"

os.system('SurveySim '+modfile+' '+templatefile+' '+obsfile+' -o '+outfile)
