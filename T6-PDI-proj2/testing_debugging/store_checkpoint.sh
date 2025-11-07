#!/bin/bash
#
number_intervals=3
time_elapsed=0

for i in $( seq 1 $number_intervals); do

sleep 6h
time_elapsed=$((time_elapsed+6))

for t in $( seq $1 $2 ); do

cd run-fssh-$t

mv "run-1.restart" trajectory_checkpoints/run-$time_elapsed-hours.restart

cd ../

done

done
