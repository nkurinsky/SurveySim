#!/bin/bash

modelfile="herschel_model.fits"
sedfile="../trunk/templates/sf_templates.fits"
obsfile="../trunk/obs/observation.fits"

echo "Updating Model File"
./update_mfile.py
echo "Running Code"
SurveySim $modelfile $sedfile $obsfile -v -o ./herschel_output.fits
