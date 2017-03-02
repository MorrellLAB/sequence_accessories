#!/bin/bash

#   Make sure we can use associative arrays
$(help declare | grep 'associative' > /dev/null 2> /dev/null) || (echo "Please use a version of BASH that supports associative arrays ('declare -A'), exiting..." >&2; exit 1)

#   Set defaults
declare -a merge_deps=(bamtools parallel) # Dependencies
OUT_DEFAULT="$(pwd -P)/merged" # Default output directory
declare -x SAMPLE_DELIMITER=','

#   Usage message
function MergeUsage() {
    echo -e "\
MergeBAM: Use BAMtools to merge BAM files together \n\
\n\
Arguments:  --sample-list=<sample_list> --name-table=<table> [--outdirectory=outdirectory] \n\
Where:      <sample_list> is a list of BAM files, named with the sample names in <table> \n\
            <table> is the sample name table (see below) \n\
            [outdirectory] is an optional output directory, defaults to ${OUT_DEFAULT} \n\
\n\
The sample name table is a whitespace-delimited table where the new sample name \n\
    is in the first column and the old sample names are in subsequent columns \n\
    Each row does not need to have the same number of columns as other rows \n\
    Lines starting with a '#' symbol are treated as comments and ignored \n\
\n\
Example: \n\
#This line is ignored \n\
MergedName1    OldName1    OldName2    OldName3 \n\
MergedName2    NameOld1    NameOld2    NameOld3    NameOld4 \n\
MergedName3    Old1        Old \n\
" >&2
    exit 1
}

#   Export the function
export -f MergeUsage

function MergeBAM() {
    local sampleList="$1" # Where is our sample list?
    local mergedname="$2" # What do we call our merged BAM file?
    local outdirectory="$3" # Where do we put our merged BAM file?
    local -a samples=($(echo $4 | tr "${SAMPLE_DELIMITER}" ' ')) # An array of old sample names
    #   Find our BAM files
    local -a bamfiles=($(grep -f <(echo "${samples[@]}" | tr ' ' '\n') "${sample_list}"))
    for s in "${bamfiles[@]}"; do [[ -f "$s" ]] || (echo "Failed to find $s" >&2; exit 1)
    #   Merge the BAM files
    (set -x; bamtools merge -list <(echo "${bamfiles[@]}" | tr ' ' '\n') > "${outdirectory}/${mergedname}_merged.bam")
}

#   Export the function
export -f MergeBAM
