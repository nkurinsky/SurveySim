#!/bin/bash

outfile="Makefile.am"

fitsfiles=$(ls *.fits)

echo "Updating Model "$outfile
echo "AUTOMAKE_OPTIONS=foreign" > $outfile
echo '' >> $outfile
echo 'modeldir = $(prefix)/@PACKAGE@/model' >> $outfile
echo "dist_model_DATA = "$fitsfiles >> $outfile

