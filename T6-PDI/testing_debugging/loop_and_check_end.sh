#!/bin/bash


#PATTERN=$1
#FILE=$2

mkdir NOT_FINISHED

for i in run-fssh-*;   #TO CHANGE!!!!
do
        echo $i
        cd $i

        if grep -q 'T I M I N G' run.log;
         then
            echo "Here are the Strings with the Pattern 'T I M I N G':"
            #echo -e "$(grep $PATTERN $FILE)\n"
            cd ../
         else
            echo "Error: The Pattern was NOT Found in run-fssh-$i"
            echo "Exiting..."
            #exit 0
            cd ../
            mv $i NOT_FINISHED
        fi
done
