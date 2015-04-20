#!/bin/bash

modelfile="jwst_model.fits"
sedfile="../trunk/templates/Kirkpatrick2015_templates.fits"

echo "Updating Model File"
./update_mfile_sim.py
echo "Running Code"
SurveySim $modelfile $sedfile -v -o ./jwst_output.fits
