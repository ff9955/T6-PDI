#!/bin/bash 

new_file="last_ener_lines.txt"

for i in {0..499}; do

cd "run-fssh-$i"
last_line=$(tail -n -1 run-1.ener)

cd ../
echo "$i : $last_line" >> $new_file

done
