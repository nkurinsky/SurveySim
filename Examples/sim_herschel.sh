#!/bin/bash

modelfile="herschel_model.fits"
sedfile="../trunk/templates/Kirkpatrick2015_templates.fits"
obsfile="../trunk/obs/observation.fits"

echo "Updating Model File"
./update_mfile.py
echo "Running Code"
SurveySim $modelfile $sedfile -v -o ./herschel_sim_output.fits
