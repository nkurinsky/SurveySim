#!/bin/bash

for field in {0..1}
do
    for lf in {0..2}
    do
	for model in {0..2}
	do
	    echo "./run_analysis.py $field $model $lf"
	    ./run_analysis.py $field $model $lf
	done
    done
done
