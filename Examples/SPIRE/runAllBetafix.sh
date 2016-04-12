#!/bin/bash

fields="0 2"
lfs="1"
betas="0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.7"

hosts=$(bhosts | grep bullet | awk '{print $1}')
echo -e Hosts:\"$hosts\"

for field in $fields
do
    for lf in $lfs
    do
	for beta in $betas
	do
	    echo "./run_analysis_betafix.py $field 2 $lf $beta"
	    bsub -W 48:00 -q medium -m "$hosts" ./run_analysis_betafix.py $field 2 $lf $beta
	done
    done
done
