#!/bin/bash

target_directory=$1

cd $target_directory

for i in $( seq 400 499 );  #TO CHANGE!!!! list of the traj number you want
do
    echo $i
    cd run-fssh-$i

    restart_file="run-1.restart.bak-1" #the name of the restart file you want to use
    updated_restart="run-1.restart.bak"

    totstep=200000 #the final timestep you want to set for the next simulation job

    step=$( grep "STEP_START_VAL" $restart_file | awk {'print $NF'} )
    states=$( grep "NUMBER_DIABATIC_STATES" $restart_file | awk {'print $NF'} )

    # Truncate data files
    cd ../
    ./../bash_scripts/truncate.sh $target_directory $i $step
    cd run-fssh-$i

    # Read diabatic coefficients at restart step number into DIAB_COEFF.include
    grep -a -A $states "i = $step" run-coeff-1.xyz | tail -$states > DIAB_COEFF.include

    # Get active state at restart step number
    active_state=$( sed -n -e "/i = $step/,/Final/ p" run-sh-1.log | tail -1 | awk {'print $4'} )

    # Remove pre-existing &WAVEFUNCTION_RESTART section from restart file
    sed '/\&WAVEFUNCTION_RESTART/,/\&END WAVEFUNCTION_RESTART/d' $restart_file > temp.txt
    mv temp.txt $restart_file

    # Define new &WAVEFUNCTION_RESTART section
    wf_section="       &WAVEFUNCTION_RESTART
         RESTART_KEY  T
         @INCLUDE   DIAB_COEFF.include
         ACTIVE_STATE_RESTART  $active_state
       &END WAVEFUNCTION_RESTART"
    
    echo "$wf_section" > SECTION.txt 

    # Insert new &WAVEFUNCTION_RESTART section into restart file
    sed '/&END AOM/r SECTION.txt' $restart_file > $updated_restart

    rm SECTION.txt

    # Update number of steps remaining
    line_to_change=$( grep " STEPS " $updated_restart ) 
    remaining_steps=$(( $totstep - $step ))

    sed -i -e "s/$line_to_change/     STEPS    $remaining_steps/g" $updated_restart

    cat run.log >> run_all.log
    mv run-1.restart.bak run.inp

    cd ../

done
