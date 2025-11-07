#!/bin/bash

atoms_per_donor=44
atoms_per_acceptor=40
first_index=284
second_index=6
total_donor_atoms=13200
total_acceptor_atoms=9600
lines_per_timestep=$(($total_donor_atoms+$total_acceptor_atoms+2))
timesteps_skipped=1
number_samples=20
line_counter=0

line_counter=$(($first_index-1))
line_counter=$((line_counter*$atoms_per_donor))
line_counter=$((3 + line_counter))

for i in $( seq 1 $number_samples ); do

current_timestep=$(($i-1))
current_timestep=$((current_timestep*$timesteps_skipped*1))

echo $(($atoms_per_donor + $atoms_per_acceptor)) >> interface-dimer-$current_timestep.txt
echo " " >> interface-dimer-$current_timestep.txt

sed -n "$line_counter, $((line_counter + $atoms_per_donor - 1))p" run-pos-1.xyz  >> interface-dimer-$current_timestep.txt

lines_skipped=$(($timesteps_skipped*$lines_per_timestep))
line_counter=$((line_counter+$lines_skipped))

done

line_counter=$(($second_index-1))
line_counter=$((line_counter*$atoms_per_acceptor))
line_counter=$(($total_donor_atoms + 3 + line_counter))

for i in $( seq 1 $number_samples ); do

current_timestep=$(($i-1))
current_timestep=$((current_timestep*$timesteps_skipped*1))

sed -n "$line_counter, $((line_counter + $atoms_per_acceptor - 1))p" run-pos-1.xyz  >> interface-dimer-$current_timestep.txt

lines_skipped=$(($timesteps_skipped*$lines_per_timestep))
line_counter=$((line_counter+$lines_skipped))

done
