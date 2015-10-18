#!/usr/bin/env python

import os
import sys
import datetime
import time

sys.path.append('/usr/local/surveysim/python')

from ModelFile import *


if(len(sys.argv) < 4):
    print "Calling Sequence: "+sys.argv[0]+" field(0=COSMOS,1=SWIRE) model(0=onlySFG,1=agn,2=composites,3=cold,4=SFG+Cold) lfForm(0=MS,1=DPL,2=S)"
    quit()
else:
    field=int(sys.argv[1])
    model=int(sys.argv[2])
    lfForm=int(sys.argv[3])

#initialize
mod=ModelFile()
mod.axis1="ColorF1F2"
mod.axis2="Flux1"

#run settings; should be same regardless of model
mod.settings['verbosity']=3
mod.annealing['temp']=.03
mod.annealing['learningRate']=0.4
mod.convergence['CI']=0.90

mod.filters[0].setID("SPIRE_250")
mod.filters[1].setID("SPIRE_350")
mod.filters[2].setID("SPIRE_500")

#set parameters for field, and input file
if(field == 0):
    simname="spire_COSMOS"
    obsfile="/usr/local/surveysim/obs/L2-COSMOS_xID250_DR2.fits"
    mod.survey['area']=4.78
    mod.filters[0].limit=8.0
    mod.filters[0].err=1.6
    mod.filters[0].compN=1.48
    mod.filters[0].compB=5.90
    mod.filters[0].compM=9.27
    mod.filters[1].limit=0.1
    mod.filters[1].err=1.32
    #mod.filters[1].compN=1.48
    #mod.filters[1].compB=4.84
    #mod.filters[1].compM=5.51
    mod.filters[2].limit=0.1
    mod.filters[2].err=1.9
elif(field == 1):
    simname="spire_Lockman-SWIRE"
    obsfile="/usr/local/surveysim/obs/L5-Lockman-SWIRE_xID250_DR2.fits"
    mod.survey['area']=15.31
    mod.filters[0].limit=9.6
    mod.filters[0].err=1.92
    mod.filters[0].compN=3.80
    mod.filters[0].compB=14.06
    mod.filters[0].compM=12.72
    mod.filters[1].limit=0.1
    mod.filters[1].err=1.58
    #mod.filters[1].compN=22.58
    #mod.filters[1].compB=13.72
    #mod.filters[1].compM=-19.83
    mod.filters[2].limit=0.1
    mod.filters[2].err=2.3
else:
    raise ValueError("Invalid field")

print mod.survey['area']

#default to all fit
mod.fitAllParams()

#remove color evolution
mod.params['cexp'].value=0
mod.params['cexp'].fixed=1
mod.params['zbc'].value=0
mod.params['zbc'].fixed=1

#set luminosity function form
if(lfForm == 0):
    mod.setLF("ModifiedSchecter")
    simname=simname+"_MS"
elif(lfForm == 1):
    mod.setLF("DoublePowerLaw")
    simname=simname+"_DPL"
elif(lfForm == 2):
    mod.setLF("Schechter")
    simname=simname+"_S"

else:
    raise ValueError("Invalid lfForm")

if(model == 0):
    simname=simname+"_onlySFG"
    fixKeys=['fa0','fcomp','fcold']
    mod.params['fa0'].value=0
    mod.params['fcomp'].value=0
    mod.params['fcold'].value=0
elif(model == 1):
    simname=simname+"_agnOnly"
    fixKeys=['fcomp','fcold']
    mod.params['fcold'].value=0
    mod.params['fcomp'].value=0
elif(model == 2):
    simname=simname+"_fComp"
    fixKeys=['fcold']
    mod.params['fcold'].value=0
    mod.params['fcomp'].value=0.5
    mod.params['fcomp'].pmin=0.01
    mod.params['fcomp'].pmax=0.99
elif(model == 3):
    simname=simname+"_fCompfCold"
    fixKeys=[]
    mod.params['fcold'].value=0.5
    mod.params['fcold'].pmin=0.01
    mod.params['fcold'].pmax=0.99
    mod.params['fcomp'].value=0.5
    mod.params['fcomp'].pmin=0.01
    mod.params['fcomp'].pmax=0.99
elif(model == 4):
    simname=simname+"_onlySFGCold"
    fixKeys=['fa0','fcomp']
    mod.params['fa0'].value=0
    mod.params['fcomp'].value=0    
else:
    raise ValueError("Invalid model")

for key in fixKeys:
    mod.params[key].fixed=1

mfile=simname+"_model.fits"
outfile=simname+"_output.fits"


#parameters below should be the same regardless of model
mod.params['Alpha'].value=3.00
mod.params['Alpha'].pmin=1.00
mod.params['Alpha'].pmax=3.50
mod.params['Alpha'].fixed=0

mod.params['Beta'].value=0.52
mod.params['Beta'].pmin=0.01
mod.params['Beta'].pmax=0.55
mod.params['Beta'].fixed=0

mod.params['Phi0'].value=-2.239
mod.params['Phi0'].pmin=-4.239
mod.params['Phi0'].pmax=-0.739
mod.params['Phi0'].fixed=0

mod.params['L0'].value=9.949
mod.params['L0'].pmin=8.449
mod.params['L0'].pmax=10.400
mod.params['L0'].fixed=0

mod.params['P'].value=-0.57
mod.params['P'].pmin=-4.57
mod.params['P'].pmax=3.43
mod.params['P'].fixed=0

mod.params['P2'].value=-2.40
mod.params['P2'].pmin=-6.40
mod.params['P2'].pmax=1.60
mod.params['P2'].fixed=0

mod.params['Q'].value=3.55
mod.params['Q'].pmin=3.05
mod.params['Q'].pmax=6.05
mod.params['Q'].fixed=0

mod.params['Q2'].value=0.80
mod.params['Q2'].pmin=-2.20
mod.params['Q2'].pmax=1.50
mod.params['Q2'].fixed=0

mod.params['zbp'].value=1.10
mod.params['zbp'].pmin=0.10
mod.params['zbp'].pmax=3.10
mod.params['zbp'].fixed=0

mod.params['zbq'].value=1.85
mod.params['zbq'].pmin=0.85
mod.params['zbq'].pmax=3.85
mod.params['zbq'].fixed=0

mod.filename=mfile
mod.run(obsfile,outfile=outfile)
