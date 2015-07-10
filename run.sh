#!/bin/bash

#a command-line means of running SurveySim -- recommended when using pre-built model files, although these can be updated by uncommenting the following lines:
#echo "Updating Model File"
#Python/update_mfile.py

sedfile="trunk/templates/Kirkpatrick2015_templates.fits"

#Comment/uncomment as needed
echo "Running SurveySim...."
#if running in fitting-mode
#modelfile="trunk/model/herschel_model.fits"
#obsfile="trunk/obs/observation.fits"
#SurveySim $modelfile $sedfile $obsfile -v -o OUTPUT/herschel_output.fits

modelfile="Examples/spire_MS_D_Fit/spire_MS_D_Fit_model.fits"
obsfile="trunk/obs/spire_cdfs-swire.fits"
outputfile="OUTPUT/spire_D_output.fits"
SurveySim $modelfile $sedfile $obsfile -v -o outputfile 

#if running in simulation-mode
#simulating Herschel
#modelfile="trunk/model/herschel_model.fits"
#SurveySim $modelfile $sedfile $obsfile -v -o OUTPUT/herschel_output.fits

#simulating JWST
#modelfile="jwst_model.fits"
#SurveySim $modelfile $sedfile -v -o OUTPUT/jwst_output.fits


