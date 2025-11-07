#!/bin/bash

path="$1/run-fssh-$2"
step=$3

#arrays of the names of files you want to truncate, if you are not printing a file that is included in this array, then remove it, otherwise it will cause an error

files=( "$path/run-coeff-1.xyz" "$path/run-pseudo-hamilt-1.xyz" "$path/run-sh-1.log" "$path/run-nace-active-1.xyz" "$path/run-adiab_pop-1.xyz" )

for file in "${files[@]}"
do
	line_number_1=$(grep -a -n "i = $step" $file | head -n 1 | cut -d: -f1) 
	if [ -z "$line_number_1" ]; then
    		echo "Error: 'i = $step' not found in file '$file'"
   		exit 1
	fi


	line_number_2=$(awk "NR > $line_number_1 { if (\$0 ~ /i = /) { print NR; exit } }" $file)
	if [ -z "$line_number_2" ]; then
    		echo "Error: 'i = $step' not found in file '$file'"
   		exit 1
	fi
	line_number_2=$((line_number_2-1))

	sed -i "${line_number_2},\$d" $file
done


files_2=( "$path/run-1.ener" )

for file in "${files_2[@]}"
do
	line_number_1=$(awk -v search_num="$step" '{if ($1 == search_num) { print NR; exit}}' $file)
	if [ -z "$line_number_1" ]; then
    		echo "Error: '$step' not found in the first column of file '$file'"
    		exit 1
	fi
	line_number_1=$((line_number_1+1))

	sed -i "${line_number_1},\$d" $file
done


