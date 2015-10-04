#!/bin/bash

#files names follow the convention:
#  - PARAMETERS file:           params_files/<testname>.txt
#  - expected OUTPUT:  output/expected/<testname>/ {files here must have}
#  - computed OUTPUT:  output/computed/<testname>/ {   the same names   }

#bin="mpirun -np 2 ./ZInCo"           # The application (from command arg)
bin="#"
diff="diff -qriad"   # Diff command, or what ever

# An array, do not have to declare it, but is supposedly faster
declare -a tests_list=("dilution" "zoom" "cascade")

# Loop the array
for testname in "${tests_list[@]}"; do
    # Padd file_base with suffixes
    file_in="params_files/$testname.txt"             # The in file
    file_expout="output/expected/$testname/"       # The out file to check against
    file_comout="output/computed/$testname/"   # The outfile from test application

    # Validate infile exists
    if [ ! -f "$file_in" ]; then
        printf "In file %s is missing\n" "$file_in"
        continue;
    fi

    printf "Testing against %s\n" "$file_in"

    # Run application, redirect in file to app, and output to out file
    "$bin $file_in > $file_in.log"

    # Execute diff
    $diff "$file_expout" "$file_comout"


    # Check exit code from previous command (ie diff)
    # We need to add this to a variable else we can't print it
    # as it will be changed by the if [
    # Iff not 0 then the files differ (at least with diff)
    e_code=$?
    if [ $e_code != 0 ]; then
            printf "TEST %s FAIL with exit code %d\n" "$testname" "$e_code"
    else
            printf "TEST %s OK!\n" "$testname"
    fi

done