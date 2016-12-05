#!/bin/env sh

#   Written by Peter L. Morrell, 18 October 2016, St. Paul, MN
#   Adapted by Paul Hoffman
#   Inspired by Tom Kono's SRA_Fetch.sh script
#   Original SRA_Fetch.sh available at:
#   https://github.com/TomJKono/Misc_Utils/blob/master/SRA_Fetch.sh

#   Set defaults
declare -a sra_deps=(lftp parallel) # Dependencies
declare -a sra_types=(experiment run sample study) # Valid SRA types
OUTPUT_DEFAULT="$(pwd -P)/SRA" # Default output directory
SRA_FTP='ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads' # Base URL for SRA

#   Usage message
function SRAUsage(){
    echo -e "\
SRADownloader: Use LFTP to download SRA files by experiment, run, sample,\n\
    or study. Optionally, use vdb-validate to validate downloaded files\n\
\n\
Arguments:  --sample-list=<sample_list> --sample-type=<sample_type>\n\
                [--outdirectory=outdirectory] [--validate]\n\
Where:      <sample_list> is a list of SRA Accessions\n\
            <sample_type> is one of 'experiment', 'run', 'sample', or 'study'\n\
            [outdirectory] is an optional output directory\n\
                (defaults to '${OUTPUT_DEFAULT}')
            [--validate] to validate downloaded SRA files with vdb-validate\n\
" >&2
    exit 1
}

#   Export the function
export -f SRAUsage

function downloadSRA() {
    local sample="$1" # What sample are we downloading?
    local ftpurl="$2" # What is the base FTP url?
    local outdirectory="$3" # Where do we store our files?
    #   Create our query
    local query="${ftpurl}/${sample:0:3}/${sample:0:6}/${sample}/"
    #   Ensure the output directory exists
    (set -x; mkdir -p "${outdirectory}")
    cd "${outdirectory}"
    #   Mirror the query directory to here
    (set -x; lftp -c "mirror ${query} $(pwd)")
    #   Make the output directory and its contents writeable
    chmod -R +w "${outdirectory}"
}

#   Export the function
export -f downloadSRA

function validateSRA() {
    local sradirectory="$1" # Where are the SRA files?
    #   Run vdb-validate
    local -a srafiles=($(find "${sradirectory}" -type f -name "*.sra"))
    echo "Validating ${#srafiles[@]} SRA files in ${sradirectory}" >&2
    parallel --verbose "vdb-validate {}" ::: "${srafiles[@]}"
}

#   Export the function
export -f validateSRA

#   Download SRA files
function SRADownloader() {
    local sampleList="$1" # Where is our sample list?
    local sampleType="$2" # What kind of SRA accession were we given?
    local outdirectory="$3" # Where do we store the output files?
    local validate="$4" # Do we validate?
    #   Assemble our base URL
    case "${sampleType}"
        experiment)
            local ftpsite="${SRA_FTP}/ByExp/sra"
            ;;
        run)
            local ftpsite="${SRA_FTP}/ByRun/sra"
            ;;
        sample)
            local ftpsite="${SRA_FTP}/BySample/sra"
            ;;
        study)
            local ftpsite="${SRA_FTP}/ByStudy/sra"
            ;;
        *)
            echo "Invalid SRA accession type: ${sampleType}" >&2; exit 1
            ;;
    esac
    #   Download files
    parallel --verbose "downloadSRA {} ${ftpsite} ${outdirectory}" :::: "${sampleList}"
    echo "Downloaded SRA files can be found at ${output}" >&2
    if [[ "${validate}" ]]; then validateSRA "${outdirectory}"; fi
}

#   Export the function
export -f SRADownloader
