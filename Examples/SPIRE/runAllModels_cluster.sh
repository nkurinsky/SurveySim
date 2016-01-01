#!/bin/bash

fields="2"
agnFracs="6.0 12.0"

hosts=$(bhosts | grep bullet | awk '{print $1}')
echo -e Hosts:\"$hosts\"

for frac in $agnFracs
do
    for field in $fields
    do
	for lf in {0..1}
	do
	    for model in {0..2}
	    do
		echo "./run_analysis_cluster.py $field $model $lf $frac"
		bsub -W 24:00 -q medium -m "$hosts" ./run_analysis_cluster.py $field $model $lf $frac
	    done
	done
    done
done
