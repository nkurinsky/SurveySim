#!/bin/bash

outfile="Makefile.am"

fitsfiles=$(ls *.fits)

echo "Updating Templates "$outfile
echo "AUTOMAKE_OPTIONS=foreign" > $outfile
echo '' >> $outfile
echo 'templatedir = $(prefix)/@PACKAGE@/templates' >> $outfile
echo "dist_template_DATA = "$fitsfiles" README" >> $outfile
