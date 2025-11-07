#!/bin/bash

NAME="run"
NCPU=1

################################################################################

CP2K_HOME="/scratch/fivanovic/flavoured-cptk-X-SH-original/cp2k"
CP2K_EXE="$CP2K_HOME/exe/local/cp2k.sopt"

cd "run-fssh-test"
  valgrind --leak-check=full --track-origins=yes --log-fd=2 $CP2K_EXE -i ${NAME}.inp > ${NAME}.log 2>&1
cd ../


################################################################################
