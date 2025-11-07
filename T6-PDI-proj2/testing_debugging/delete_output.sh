#!/bin/bash
#
for i in {1..399}; do

        cd "run-fssh-$i"
        rm CT_DISTANCES.include
        cd ../

done
