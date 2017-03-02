#!/bin/bash

#   Set some defaults
declare -a cov_deps=(samtools awk parallel)
OUTPUT_DEFAULT="$(pwd -P)/Coverage"
PROJECT_DEFAULT='SimpleCoverage'

#   A usage message
function covUsage() {
    echo -e "\
SimpleCoverage: Use SAMTools to calculate coverage of a BAM file.\n\
\n\
Arguments:  --sample-list=<sample_list> [--genome-size=<genome_size>] [--project=<project>] [--outdirectory=<outdirectory>] \n\
Where:      <sample_list>   is a list of indexed BAM files \n\
            [genome_size]   is an optional genome size will calculate if not specified) \n\
            [outdirectory]  is an optional output directory, defaults to ${OUT_DEFAULT} \n\
            [project]       is an optional name for the output file, defaults to '${PROJECT_DEFAULT}' \n\
" >&2
    exit 1
}

#   Export the function
export -f covUsage

#   A function to validate integers
function isInteger() {
    local myNumber="$1"
    integer='^[0-9]+$'
    [[ "$myNumber" =~ "$re" ]] || (echo "${myNumber} is not an integer, please specify an integer" >&2; exit 1)
}

#   Export the function
export -f isInteger

#   A function to calculate genome size
function genomeSize() {
    local bamfile="$1" # What BAM file are we working with?
    local size="$(samtools view -H ${bamfile} \
        | grep -oE "LN:[[:digit:]]+" \
        | cut -f 2 -d ':' \
        | awk '{sum+=$1} END {print sum}')"
    echo "${size}"
}

#   Export the function
export -f genomeSize

#   A function to get coverage
function getCoverage() {
    local bamfile="$1"
    local genomesize="${2:-$(genomeSize ${bamfile})}"
    isInteger "${genomesize}"
    local bamname="$(basename ${bamfile} .bam)"
    #   Get the number of reads
    local numreads="$(samtools depth ${bamfile} \
        | awk '{sum+=$3} END {print sum}')"
    #   Get the raw coverage
    local rawcov="$(echo ${numreads} / ${genomesize} | bc -l)"
    #   Format it to look pretty
    local coverage="$(printf '%.3f' ${rawcov})"
    echo -e "${bamname}\t${coverage}"
}

#   Export the function
export -f getCoverage

#   A way to simply calculate coverage
function SimpleCoverage() {
    local bamlist="$1"
    local outdirectory="$2"
    local project="$3"
    local target="$4"
    parallel --verbose "getCoverage {} ${target}" > "${outdirectory}/${project}.coverage" :::: "${bamlist}"
}

#   Export the function
export -f SimpleCoverage
