#!/usr/bin/env bash
# https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
set -o pipefail
# Validate a FASTQ file. Remove lines that are invalid.

# Parse user provided options and their arguments.
# https://gist.github.com/c-garcia/95e488e974f207f3afa95fca2fdf683b
# https://www.redhat.com/sysadmin/arguments-options-bash-scripts
while getopts ":i:o:" opt; do
case $opt in
    i) fastq_to_check=$OPTARG;;
    o) polished_fastq=$OPTARG;;
    \?) # Invalid option
        echo "Error: Invalid option"
        exit;;
esac
done

# File checks
# Check that fastq_to_check exists.
# if ! [ -f "${fastq_to_check}" ]
# then
#     printf "${fastq_to_check} cannot be found.\n"
#     exit 1
# fi
if [ -f "${polished_fastq}.log" ]
then
    # ${polished_fastq}.log exists
    rm ${polished_fastq}.log
fi

#####################################################################
# Check for specific file format problems.
# Note that fastq_info will terminate upon finding the first error.
# Using the trick here,
# https://stackoverflow.com/a/2342841/8423001
# get the line with the problem.
 
# Work with a temporary file.  Iteratively delete problematic lines
# from the temporary file until the temporary file no longer has 
# any problematic lines.

cp ${fastq_to_check} ${polished_fastq}.tmp
# last_lines_in_problem_reads will usually give the last line in the 4-line 
# read sequence with the problem.
last_line_with_new_problem=$(fastq_info ${polished_fastq}.tmp 2>&1 1> ${polished_fastq}.log \
    | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+" -)

while [[ -n "${last_line_with_new_problem}" ]]
do
    # last_line_with_new_problem is of non-zero length, i.e.
    # we found a problem for a read.  
    # Delete the problematic line.
    # https://stackoverflow.com/a/48744059/8423001
    # https://stackoverflow.com/a/26569006/8423001
    # https://stackoverflow.com/a/15978536/8423001
    # https://stackoverflow.com/a/26727351/8423001
    sed "${last_line_with_new_problem}d;" "${polished_fastq}.tmp" \
        > ${polished_fastq}.pre_tmp && mv ${polished_fastq}.pre_tmp ${polished_fastq}.tmp

    # broken pipe problems: https://superuser.com/a/642932/1774660
    last_line_with_new_problem=$(fastq_info ${polished_fastq}.tmp 2>&1 1>> ${polished_fastq}.log \
        | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+" -)
done

mv ${polished_fastq}.tmp ${polished_fastq}