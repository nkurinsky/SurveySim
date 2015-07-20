#!/bin/bash

outfile="Makefile.am"

fitsfiles=$(ls *.fits)

echo "Updating Obs "$outfile
echo "AUTOMAKE_OPTIONS=foreign" > $outfile
echo '' >> $outfile
echo 'obsdir = $(prefix)/@PACKAGE@/obs' >> $outfile
echo "dist_obs_DATA = "$fitsfiles" OBSINFO" >> $outfile
