#!/bin/bash

declare -a sub_deps(seqtk parallel)

function subsampleUsage() {
    echo -e "\
" >&2
    exit 1
}

export -f subsampleUsage

function subsample() {
    local sample="$1" # What sample are we working with?
    local seed="$2" # What is our seed?
    local fraction="$3" # What fraction of reads are we keeping?
    local outdirectory="$4" # Where do we put our finished file?
    #   Ensure we have an output directory
    mkdir -p "${outdirectory}"
    #   Figure out extensions
    local extension="$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" # Get the last bit of the extension
    if [[ "${extension}" == 'gz' || "${extension}" == 'bz2' ]] # If we have a gzipped or bzipped2 FASTQ file
    then
        local fullextension="$(echo ${sample} | rev | cut -f 1,2 -d '.' | rev)" # Then get the last two extensions
    else
        local fullextension="${extension}" # Otherwise, just the last
    fi
    #   Figure out our sample name
    local samplename="$(basename ${sample} .${fullextension})" # Remove the directory and full extension from our sample name
    local decomp="${outdirectory}/${samplename}_PIPE" # Create a name for decompression (only used for pipes)
    rm -f "${decomp}" # Remove any possible pipes
    mkfifo "${decomp}" # Make a pipe for it
    #   Create an output name
    local outbase="${outdirectory}/${samplename}_${fraction}.fastq" # Base output name
    if [[ "${extension}" == 'gz' ]] # If we started with a gzipped file
    then
        local outname="${outbase}.gz" # Create a gzipped file
        local compression='gzip' # We recompress with gzip
        gzip -cd "${sample}" > "${decomp}" & # Decompress to the pipe
    elif [[ "${extension}" == 'bz2' ]] # If we stared with a bzipped2 file
    then
        local outname="${outbase}.bz2" # Create a bzipped2 file
        local compression='bzip2' # Recompress with bzip2
        bzip2 -cd "${sample}" > "${decomp}" & # Decompress to the pipe
    else
        local outname="${outbase}" # Output name is same as base output name
        local compression='tee /dev/null' # Don't recompress
        cat "${sample}" > "${decomp}" & # Write sample to pipe for easier cleanup later
    fi
    (set -x; seqtk sample -s "${seed}" "${decomp}" "${fraction}" | "${compression}" > "${outname}")
    find "${outdirectory}" -type p -name "$(basename ${decomp})" -exec rm -f {} \; # Remove the pipe (if used)
}

function SubsampleFastq() {
    local sampleList=$1

}

export -f SubsampleFastq
