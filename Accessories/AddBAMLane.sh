#!/bin/bash

set -eo pipefail

#   Some defaults
declare -a lane_deps=(samtools parallel)
OUTPUT_DEFAULT="$(pwd -P)/withLane"

#   Usage message
function LaneUsage() {
    echo -e "\
" >&2
    exit 1
}

#   Export the function
export -f LaneUsage


#   Read group IDs cannot have white space
function findReadGroups() {
    local bamfile="$1"
    declare -a readGroups=($(samtools view "${bamfile}" | grep -oE 'RG(:Z)?:\S+' | cut -f 3 -d ':' | sort -u))
    echo ${readGroups[@]}
}

function parseRGHeader() {
    local bamfile="$1"
    declare -a RGIDs=($(samtools view -H "${bamfile}" | grep -oE '@RG.*ID:\S+' | cut -f 2 -d ':'))
}
