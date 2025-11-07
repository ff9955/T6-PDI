#!/bin/bash

atoms_per_donor=44
atoms_per_acceptor=40
total_donor_molecules=300
total_donor_atoms=$((atoms_per_donor*$total_donor_molecules))

donor_indices=()
acceptor_indices=()
atom_counter=0

test_geom_file="defect_production_run.xyz"
file_to_read="final_nvt_geom_defect.xyz"

for j in {1..300..30}; do

initial_index=$((j+8))
final_index=$((initial_index+10))

for k in $( seq $initial_index $final_index ); do
donor_indices+=( $k )
done  

done
echo ${donor_indices[@]}

for j in {301..540..12}; do

initial_index=$((j+3))
final_index=$((initial_index+4))

for k in $( seq $initial_index $final_index ); do
acceptor_indices+=( $k )
done

done
echo ${acceptor_indices[@]}

donor_indices=( 14 15 44 45 74 75 104 105 134 135 164 165 194 195 224 225 254 255 284 285 )
acceptor_indices=( 306 329 341 353 365 377 389 401 413 425 437 449 461 473 485 497 509 521 533 )

for element in ${donor_indices[@]}; do

line_counter=3

index_multiplier=$((element-1))
lines_skipped=$((index_multiplier*$atoms_per_donor))
line_counter=$((line_counter + $lines_skipped))

sed -n "$line_counter, $((line_counter + $atoms_per_donor - 1))p" $file_to_read >> $test_geom_file

atom_counter=$((atom_counter + atoms_per_donor))

done

for element in ${acceptor_indices[@]}; do

line_counter=$((total_donor_atoms+3))

index_multiplier=$((element-1-$total_donor_molecules))
lines_skipped=$((index_multiplier*$atoms_per_acceptor))
line_counter=$((line_counter + $lines_skipped))

sed -n "$line_counter, $((line_counter + $atoms_per_acceptor - 1))p" $file_to_read >> $test_geom_file

atom_counter=$((atom_counter + atoms_per_acceptor))

done

dummy_file="dummy_file.xyz"
echo $atom_counter >> $dummy_file
echo " " >> $dummy_file

cat $test_geom_file >> $dummy_file
rm $test_geom_file
mv $dummy_file $test_geom_file
