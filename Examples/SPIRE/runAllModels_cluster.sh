#!/bin/bash

fields="0 2"
lfs="0 1"
agnFracs="6.0 12.0"

hosts=$(bhosts | grep bullet | awk '{print $1}')
echo -e Hosts:\"$hosts\"

for field in $fields
do
    for lf in $lfs
    do
	echo "./run_analysis_cluster.py $field 2 $lf"
	#bsub -W 24:00 -q medium -m "$hosts" ./run_analysis_cluster.py $field 2 $lf
    done
done
