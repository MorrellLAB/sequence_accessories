#!/bin/env bash

#   Written by Peter L. Morrell, 26 July 2016, St. Paul, MN
#   Updated 27 October 2016, 29 November 2016
#   Modified by Paul Hoffman

#   Dependencies for DumpFastq
declare -a dump_deps=(parallel fastq-dump)
OUTPUT_DEFAULT="$(pwd -P)/FASTQ"

#   Create a usage message
function dumpUsage() {
    echo -e "\
$(basename $0) DumpFastq: Use fastq-dump from the SRA Toolkit to dump \n\
        SRA archive files (.sra) to gzipped FASTQ files (.fastq.gz) \n\
\n\
Usage:  $(basename $0) DumpFastq --sample-list=sample_list [--outdirectory=outdirectory] [--paired] \n\
Where:  <sample_list> is a list of SRA files to be dumped to FASTQ format \n\
        [outdirectory] is an optional output directory (defaults to '${OUTPUT_DEFAULT}') \n\
        [--paired] is passed when SRA files should be dumped to paired-end FASTQ files \n\
" >&2
    exit 1
}

#   Export the function
export -f dumpUsage

#   DumpFastq
function DumpFastq() {
    local sampleList="$1" # Where are our samples?
    local outdir="$2" # Where do we put our outputs?
    local paired="$3" # Are our samples paired end?
    #   Make an out directory
    (set -x; mkdir -p "${outdir}")
    #   Set fastq-dump options
    local dumpOpts="-I -F --gzip --outdir ${outdir}"
    if ${paired}; then local dumpOpts="--split-files ${dumpOpts}"; fi
    #   Run fastq-dump
    parallel --verbose "fastq-dump ${dumpOpts} {}" :::: "${sampleList}"
    echo "Gzipped FASTQ files can be found at ${outdir}" >&2
}

#   Export the function
export -f DumpFastq
