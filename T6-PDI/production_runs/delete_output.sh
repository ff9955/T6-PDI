#!/bin/bash
#
for i in {0..499}; do

        cd "run-fssh-$i"
        rm *run-*
        cd ../

done
