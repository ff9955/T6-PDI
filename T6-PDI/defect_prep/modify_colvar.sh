#!/bin/bash

colvar_counter=0

#Check if a file was provided
if [ -z "$1" ]; then
	  echo "Usage: $0 <filename>"
	    exit 1
fi

input_file="$1"
output_file="defect_COLVAR.include"

# Clear or create the output file
 > "$output_file"

 # Loop through each line in the input file
 while IFS= read -r line; do
   # Process each word separately
     modified_line=""
       for word in $line; do
           if [[ "$word" =~ ^-?[0-9]+$ ]]; then
                 # If it's an integer, subtract 40
		 
		 if [ $word -eq 1 ]; then
			 colvar_counter=$(( colvar_counter + 1 ))
		 elif [ $word -lt  13883 ]; then
			 :
		 else
			 word=$(( word - 40 ))
		 fi

           fi
                               # Append to the modified line
           modified_line="$modified_line $word"
       done

       # Write the modified line to the output file (trim leading spaces)
         echo "${modified_line## }" >> "$output_file"
         done < "$input_file"

         echo "Subtraction completed. Check $output_file for the results."
	 echo "colvar counter:" $colvar_counter
