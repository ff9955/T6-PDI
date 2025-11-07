#!/bin/bash

line=1
step=0
total_lines=$(wc -l 284-306.xyz | cut -d ' ' -f1)

while [ $line -lt $total_lines ]
do

sed -n "$line, $((line+85))p" 284-306.xyz >> PDI-T6-dimer-$step.txt

line=$((line+86))
step=$((step+1))

done
