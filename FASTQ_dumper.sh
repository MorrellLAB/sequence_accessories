#!/bin/env bash

#PBS -l mem=1000mb,nodes=1:ppn=16,walltime=72:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q mesabi

set -euo pipefail

#    Peter L. Morrell, 26 July 2016, St. Paul, MN
#    Updated 27 October 2016, 29 November 2016
#    Dependencies: SRA Toolkit

module load sratoolkit
module load parallel

#    directory for output of fastq.gz files
WORKING=/panfs/roc/scratch/pmorrell/MBE_Barley

#   initalize the array that will hold a list of SRA files
declare -a SRA=($(find ${WORKING} -maxdepth 1 -name '*.sra' -type f))
echo "${#SRA[@]} samples"

parallel --verbose "fastq-dump --split-files -I -F --gzip {} --outdir ${WORKING}" ::: ${SRA[@]}

