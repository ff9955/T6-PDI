#!/bin/bash

cd ../../optimise_timestep_for2D/0.03fs

for j in {0..99}; do

cd run-fssh-$j

grep "CASE" run_all.log >> cases.log

grep "fs" run_all.log >> all_case_times.log
grep "=" all_case_times.log >> case_times.log
rm all_case_times.log

rm run.log

cd ../

done
