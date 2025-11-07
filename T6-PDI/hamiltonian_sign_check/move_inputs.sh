#!/bin/bash
#
for i in $( seq 0 9); do

        mkdir step-$i
        mv *$i.* step-$i
done

