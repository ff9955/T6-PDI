#!/bin/bash
#
counter=0

for j in $(seq 0 9); do

cd step-$j

sbatch submit-short-$j.sh
counter=$((counter+1))

cd ../

if [ $counter -eq 4 ]; then
	counter=0
	sleep 20m
fi

done
