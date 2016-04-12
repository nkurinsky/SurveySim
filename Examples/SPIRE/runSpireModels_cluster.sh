#!/bin/bash

fields="0"
lfs="1"

hosts=$(bhosts | grep bullet | awk '{print $1}')
echo -e Hosts:\"$hosts\"

for field in $fields
do
    for lf in $lfs
    do
	#echo "./run_analysis_cluster.py $field 2 $lf"
	#bsub -W 48:00 -q medium -m "$hosts" ./run_analysis_cluster.py $field 2 $lf
	echo "./run_analysis_cluster_fixz.py $field 2 $lf"
        bsub -W 48:00 -q medium -m "$hosts" ./run_analysis_cluster_fixz.py $field 2 $lf
    done
done
