#!/bin/bash

#   Set some defaults
declare -a list_deps=(find)
OUT_DEFAULT="$(pwd -P)/samples.txt"

#   Usage message
function ListUsage() {
    echo -e "\
ListGenerator: Create a sample list for use with sequence_handling and sequence_accessories \n\
\n\
Arguments:  --sample-directory=<directory> [--extension=extension] [--output-list=output_list] [--append-list=append_list] [--follow-tree] \n\
Where:      <directory> is the directory where your samples are \n\
            [extension] is an optional extension to limit the search to, defaults to all files \n\
            [output_list] is an optional name to create a new list, defaults to ${OUT_DEFAULT} \n\
                Note: \`--output-list' is mutually exclusive with \`--append-list' \n\
            [append_list] is an optional preexisting sample list to append to \n\
            [--follow-tree] to look at subdirectories within <directory> \n\
                Otherwise, ListGenerator will look in <directory> and NO subdirectories \n\
" >&2
    exit 1
}

#   Export the function
export -f ListUsage
