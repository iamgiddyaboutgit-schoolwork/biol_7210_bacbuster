#!/bin/bash
fastq_to_check="../../testing_data/sequencing_reads/test_fastq_3.fq"
num_lines_in_file=$(wc -l ${fastq_to_check} | cut -f 1 -d " ")

# Note that fastq_info will terminate upon finding the first error.
# Using the trick here,
# https://stackoverflow.com/a/2342841/8423001
# get the line with the problem.

# last_line_for_problem gives the last line in the 4-line sequence with the problem.
# Problems that we use a regex to find are "duplicated sequence" and "sequence and 
# quality score lengths don't match".
# We need to delete all 4 lines with the problem in the .fq file.
# The [0] allows us to assign to a slot in an array.
# https://www.gnu.org/software/bash/manual/html_node/Arrays.html
last_line_for_problem[0]=$(fastq_info ${fastq_to_check} 2>&1 > /dev/null \
    | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)

# https://stackoverflow.com/a/14900557/8423001
# https://stackoverflow.com/q/3869072/8423001
# Check if last_line_for_problem is non-empty.
if [[ -n "${last_line_for_problem[0]}" ]] && (( "${last_line_for_problem[0]}" < "${num_lines_in_file}" )); then
    # https://unix.stackexchange.com/a/306141
    line_after_problem=$((${last_line_for_problem}+1))

    # Splice out problematic lines from fastq_to_check.
    # The + sign in the -n option to tail tells tail to output
    # starting at the given line number.
    last_line_for_problem+=$(tail -n +${line_after_problem} ${fastq_to_check} \
        | fastq_info - 2>&1 > /dev/null \
        | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)
fi
my_var=2
my_var=$((${my_var}+1))
echo ${my_var}
echo ${line_before_problem}
echo $(head -n ${line_before_problem} ${fastq_to_check} \
        | fastq_info - )
# # Now, iterate, taking lines out of polished.fq as necessary.
# # More information on iterating through an array: https://www.geeksforgeeks.org/bash-scripting-array/
# last_line_for_problem_again=$(fastq_info polished.fq 2>&1 > /dev/null \
#     | grep -P --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)
# while [[ -n "$last_line_for_problem_again" ]]
# do
#     line_before_problem_again=$((${last_line_for_problem_again}-4))
#     line_after_problem_again=$((${last_line_for_problem_again}+1))
#     head -n ${line_before_problem_again} polished.fq > "polished.fq"
#     tail -n +${line_after_problem_again} polished.fq >> "polished.fq"
# done


