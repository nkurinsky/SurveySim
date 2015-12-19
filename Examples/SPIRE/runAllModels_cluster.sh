#!/bin/bash

fields="2"

for field in $fields
do
    for lf in {0..1}
    do
	for model in {0..2}
	do
	    echo "./run_analysis_cluster.py $field $model $lf"
	    bsub -W 960 ./run_analysis_cluster.py $field $model $lf
	done
    done
done
