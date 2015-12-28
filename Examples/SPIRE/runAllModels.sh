#!/bin/bash

export PATH=$PATH:/cluster/tufts/sajina_lab/surveysim/trunk/src

sedfile=/cluster/tufts/sajina_lab/surveysim/trunk/templates/default_templates.fits
obsfile=/cluster/tufts/sajina_lab/surveysim/trunk/obs/L2-COSMOS_xID24_DR3.fits

#mips runs
SurveySim A_mips_model.fits $sedfile $obsfile -o A_mips_output.fits
SurveySim B_mips_model.fits $sedfile $obsfile -o B_mips_output.fits
SurveySim C_mips_model.fits $sedfile $obsfile -o C_mips_output.fits
SurveySim D_mips_model.fits $sedfile $obsfile -o D_mips_output.fits
SurveySim E_mips_model.fits $sedfile $obsfile -o E_mips_output.fits
SurveySim F_mips_model.fits $sedfile $obsfile -o F_mips_output.fits

#Spire runs
obsfile=/cluster/tufts/sajina_lab/surveysim/trunk/obs/L2-COSMOS_xID250_DR2.fits

SurveySim A_spire_model.fits $sedfile $obsfile -o A_spire_output.fits
SurveySim B_spire_model.fits $sedfile $obsfile -o B_spire_output.fits
SurveySim C_spire_model.fits $sedfile $obsfile -o C_spire_output.fits
SurveySim D_spire_model.fits $sedfile $obsfile -o D_spire_output.fits
SurveySim E_spire_model.fits $sedfile $obsfile -o E_spire_output.fits
SurveySim F_spire_model.fits $sedfile $obsfile -o F_spire_output.fits
