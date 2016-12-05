#!/bin/env sh

#   Written by Peter L. Morrell, 18 October 2016, St. Paul, MN
#   Adapted by Paul Hoffman
#   Inspired by Tom Kono's SRA_Fetch.sh script
#   Original SRA_Fetch.sh available at:
#   https://github.com/TomJKono/Misc_Utils/blob/master/SRA_Fetch.sh

#    Path to file listing SRA files to download
SRA_FILES=${HOME}/scratch/IPK_Landraces/IPK_landraces_SRA.txt

#   Directory where SRA files will be downloaded
OUTPUT=${HOME}/scratch/IPK_Landraces

#   Make sure the file exists
if ! [[ -f "${SRA_FILES}" ]]
    then echo "Failed to find ${SRA_FILES}, exiting..." >&2
    exit 1
    fi

#   Make the array using command substitution
declare -a SRA_ARRAY=($(cat "${SRA_FILES}"))

#   Print the values of the array to screen
printf '%s\n' "${SRA_ARRAY[@]}"

#   location of SRA download script, from Tom Kono's Misc_Utils
SRA_FETCH=${HOME}/Apps/TKono/Misc_Utils/SRA_Fetch.sh

#   iterate over every each of the run numbers in a lit of SRA files
#   and download to specified directory
#   in SRA_Fetch -r run #, -e experiment #, -p sample #, -s study #
#   -d option specifies output directory
for i in "${SRA_ARRAY[@]}"
    do bash $SRA_FETCH -r $i -d $OUTPUT
    done

declare -a sra_deps=()

function SRAUsage(){
    echo -e "\
" >&2
    exit 1
}

export -f SRAUsage





