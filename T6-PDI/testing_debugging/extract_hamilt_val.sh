#!/bin/bash

first_element=10
second_element=11
coupling_array=()
current_step=0
relevant_step=1
steps_skipped=1

while IFS="   " read -r part1 part2 part3; do


if [ "$part1" == "$first_element" ]; then

	if [ "$part2" == "$second_element" ]; then

		current_step=$((current_step+1))

		if [ $current_step -eq $relevant_step ]; then
		
		relevant_step=$(($steps_skipped+relevant_step))
		coupling_array+=("$part3")

		fi

	fi

fi

done < run-pseudo-hamilt-1.xyz

for element in ${coupling_array[@]}; do
  echo $element >> PDI_XT_couplings.xyz 
done
