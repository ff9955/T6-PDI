#!/bin/bash

atoms_per_donor=44
atoms_per_acceptor=40
total_donor_molecules=300
total_donor_atoms=$((atoms_per_donor*$total_donor_molecules))

donor_indices=()
acceptor_indices=()
atom_counter=0

test_geom_file="final-nvt-subset-defect.xyz"
file_to_read="final_nvt_geom_defect.xyz"

for j in {1..300..30}; do

initial_index=$((j+8))
final_index=$((initial_index+10))

for k in $( seq $initial_index $final_index ); do
donor_indices+=( $k )
done  

done
echo ${donor_indices[@]}

#first PDI column
starting_pdi=301
initial_index=$((starting_pdi+3))
final_index=$((initial_index+4))
for k in $( seq $initial_index $final_index ); do
acceptor_indices+=( $k )
done

#2nd PDI column where there's one missing
starting_pdi=313
initial_index=$((starting_pdi+3))
final_index=$((initial_index+3))
for k in $( seq $initial_index $final_index ); do
acceptor_indices+=( $k )
done


for j in {324..539..12}; do

initial_index=$((j+3))
final_index=$((initial_index+4))

for k in $( seq $initial_index $final_index ); do
acceptor_indices+=( $k )
done

done
echo ${acceptor_indices[@]}

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
