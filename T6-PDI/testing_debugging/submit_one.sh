#!/bin/bash

NAME="run"
NCPU=1

################################################################################

CP2K_HOME="/scratch/fivanovic/flavoured-cptk-X-SH_uncommented/flavoured-cptk-X-SH/cp2k"
CP2K_EXE="$CP2K_HOME/exe/local/cp2k.sopt"

cd "run-fssh-test"
  $CP2K_EXE -i ${NAME}.inp > ${NAME}.log &
cd ../


################################################################################
