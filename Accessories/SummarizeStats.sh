#!/bin/bash

#   Dependencies for SummarizeStats
declare -a idx_deps=(datamash parallel samtools)
PROJECT_DEFAULT='STATS'

#   Create a usage message
function idxUsage() {
    echo -e "\
$(basename $0) SummarizeStats: Use SAMTools idxstats to generate simple coverage \n\
    statistics and create a summary statiscs report \n\
\n\
Usage:  $(basename $0) SummarizeStats --sample-list=<sample_list> [--project=project] \n\
Where:  <sample_list>   is a list of indexed BAM files\n\
        [project]       is an optional name for the output file, defaults to '${PROJECT_DEFAULT}' \n\
" >&2
    exit 1
}

#   Export the function
export -f idxUsage

#   A function to check for and index BAM files
function checkIndex() {
    local bamfile="$1" # Where is our BAM file?
    if ! [[ -f "${bamfile}" ]]; then echo "Cannot find BAM file, exiting..." >&2; exit 1; fi
    #   Check for an index file
    if ! [[ -f "${bamfile/bam/csi}" \
        || -f "${bamfile/bam/bai}" \
        || -f "${bamfile}.csi" \
        || -f "${bamfile}.bai" ]]
    then
        #   If no index, then generate a CSI index
        echo "Indexing ${bamfile}" >&2
        samtools index -c "${bamfile}"
    fi
}

#   Export the function
export -f checkIndex

#   A function to sum the sequence length, number of mapped reads, and number of unmapped reads for a single sample
function getCounts() {
    local sample="$1" # What sample are we working on?
    local sampleExtension=".$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" # Get the extension
    local sampleName="$(basename ${sample} ${sampleExtension})" # Remove the extension and file path
    local sequenceLength=$(cat ${sample} | head -$(($(wc -l < ${sample})-1)) | datamash sum 2) # Get the total sequence length
    local mappedReads=$(cat ${sample} | head -$(($(wc -l < ${sample})-1)) | datamash sum 3) # Get the total number of mapped reads
    local unmappedReads=$(cat ${sample} | head -$(($(wc -l < ${sample})-1)) | datamash sum 4) # Get the total number of unmapped reads
    echo -e "${sampleName}\t${sequenceLength}\t${mappedReads}\t${unmappedReads}" # Return the values, tab-delimeted
}

#   Export the function
export -f getCounts

#   Run SummarizeStats
function SummarizeStats() {
    local sampleList="$1" # List of BAM files
    local project="$2" # What do we call our output?
    #   Set output options
    local stats_dir="$(dirname ${sampleList})/idxstats" # Directory for statistics
    local output="${stats_dir}/${project}_allStats.txt" # Summarized output file
    #   Make the output directory
    mkdir -p "${stats_dir}"
    #   Check for and index BAM files
    echo "Checking indexes for BAM files..." >&2
    parallel --verbose checkIndex :::: "${sampleList}"
    #   Generate statisitcs for BAM files
    echo "Generating statistics for BAM files..." >&2
    parallel --verbose --xapply "samtools idxstats {1} > ${stats_dir}/{1/.}.idxstats" :::: "${sampleList}"
    #   Summarize the statiscs
    local -a idxstats=($(find "${stats_dir}" -maxdepth 1 -type f -name "*.idxstats")) # Find the idxstats files
    echo "Summarizing statistics..." >&2
    parallel getCounts {} > "${output}" ::: "${idxstats[@]}"
}

#   Export the function
export -f SummarizeStats
