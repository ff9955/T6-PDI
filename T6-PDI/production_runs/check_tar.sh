#!/bin/bash
tar_output="../physopt_0.01fs/tar-$1-$2.log"

for j in $( seq $1 $2 ); do

	trajectory_info=$(grep "run-fssh-$j" $tar_output | wc -l)

	if [ $trajectory_info == "0" ]; then
		echo "$j missing"
	fi

done
