#!/bin/bash

for field in "0 2"
do
    for lf in "0 1"
    do
	for model in {0..2}
	do
	    echo "./run_analysis_cluster.py $field $model $lf"
	    bsub -q long ./run_analysis_cluster.py $field $model $lf
	done
    done
done
