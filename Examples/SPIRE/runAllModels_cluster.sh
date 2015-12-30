#!/bin/bash

fields="2"

hosts=$(bhosts | grep bullet | awk '{print $1}')
echo -e Hosts:\"$hosts\"

for field in $fields
do
    for lf in {0..1}
    do
	for model in {0..2}
	do
	    echo "./run_analysis_cluster.py $field $model $lf"
	    bsub -W 24:00 -q medium -m "$hosts" ./run_analysis_cluster.py $field $model $lf
	done
    done
done
