#!/usr/bin/env bash
fastq_to_check="../../testing_data/sequencing_reads/test_fastq_3.fq"
num_lines_in_file=$(wc -l ${fastq_to_check} | cut -f 1 -d " ")

#####################################################################
# Check for any problems.
# problem_counter=0
# last_lines_in_problem_reads[${problem_counter}]=$(fastq_info ${fastq_to_check} 2>&1 > /dev/null \
#     | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)
# echo ${last_lines_in_problem_reads}
# Iteratively check for additional problems if there was a 
# previous problem.


# while (( "${line_after_problem}" < "${num_lines_in_file}" ))
# do
#     last_lines_in_problem_reads[${problem_counter}]=$(tail -n +${line_after_problem} ${fastq_to_check} \
#         | fastq_info - 2>&1 > /dev/null \
#         | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)

#     line_after_problem=$((${last_lines_in_problem_reads[${problem_counter}]}+1))

#     ((problem_counter+=1))
#     ##### loop check ####
#     echo "line_after_problem: ${line_after_problem}"
#     echo "num_lines_in_file: ${num_lines_in_file}"
# done

# Note that fastq_info will terminate upon finding the first error.
# Using the trick here,
# https://stackoverflow.com/a/2342841/8423001
# get the line with the problem.

# last_lines_in_problem_reads gives the last line in the 4-line sequence with the problem.
# Problems that we use a regex to find are "duplicated sequence" and "sequence and 
# quality score lengths don't match".
# We need to delete all 4 lines with the problem in the .fq file.
# The [0] allows us to assign to a slot in an array.
# https://www.gnu.org/software/bash/manual/html_node/Arrays.html
# last_lines_in_problem_reads[0]=$(fastq_info ${fastq_to_check} 2>&1 > /dev/null \
#     | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)

# # https://stackoverflow.com/a/14900557/8423001
# # https://stackoverflow.com/q/3869072/8423001
# # Check if last_lines_in_problem_reads is non-empty and last_lines_in_problem_reads[0] < num_lines_in_file.
# if [[ -n "${last_lines_in_problem_reads[0]}" ]] && (( "${last_lines_in_problem_reads[0]}" < "${num_lines_in_file}" )); then
#     # https://unix.stackexchange.com/a/306141
#     line_after_problem=$((${last_lines_in_problem_reads}+1))

#     # Splice out problematic lines from fastq_to_check.
#     # The + sign in the -n option to tail tells tail to output
#     # starting at the given line number.
#     next_last_line=$(tail -n +${line_after_problem} ${fastq_to_check} \
#         | fastq_info - 2>&1 > /dev/null \
#         | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)

#     # Was there another problem (of the kind that we are looking for) in the file?
#     if [[ -n "${next_last_line}" ]]; then
#         last_lines_in_problem_reads+=$(${next_last_line})
#     fi
# fi






# # Now, iterate, taking lines out of polished.fq as necessary.
# # More information on iterating through an array: https://www.geeksforgeeks.org/bash-scripting-array/
# last_lines_in_problem_reads_again=$(fastq_info polished.fq 2>&1 > /dev/null \
#     | grep -P --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)
# while [[ -n "$last_lines_in_problem_reads_again" ]]
# do
#     line_before_problem_again=$((${last_lines_in_problem_reads_again}-4))
#     line_after_problem_again=$((${last_lines_in_problem_reads_again}+1))
#     head -n ${line_before_problem_again} polished.fq > "polished.fq"
#     tail -n +${line_after_problem_again} polished.fq >> "polished.fq"
# done

# A loop can be used to get a running sum of the elements of 
# last_lines_in_problem_reads.  From this, we can infer
# the line numbers in the original file.
# The idea is that if we want to do all of our deleting of lines
# in the file in memory, then we need to know the actual line
# numbers that need to be deleted ahead of time.  
# We can easily get an array of line numbers using fastq_info
# that are based on lines since the last recursive call.
# To change this array to absolute line numbers, we just
# need to perform some simple arithmetic.
# Declare an array.
declare -a last_lines_in_problem_reads_reformatted 
# TODO: delete
declare -a last_lines_in_problem_reads

last_lines_in_problem_reads=(4 8 16 24) # TODO:for testing only
num_last_lines_in_problem_reads=${#last_lines_in_problem_reads[@]}

if  [[ -n "${last_lines_in_problem_reads[0]}" ]]; then
    # last_lines_in_problem_reads[0] is of non-zero length
    last_lines_in_problem_reads_reformatted[0]=${last_lines_in_problem_reads[0]}
fi

j=1
# After this loop, last_lines_in_problem_reads_reformatted will hold
# line numbers positioned within the original file.
for (( i=1 ; i<${num_last_lines_in_problem_reads} ; i++ )); do
    last_lines_in_problem_reads_reformatted[${j}]=$((${last_lines_in_problem_reads_reformatted[$((${j} - 1))]} + ${last_lines_in_problem_reads[i]}))
    j=$((${j}+1))
    
done

declare -a all_lines_to_delete
index_for_all_lines_to_delete=0
for line_num in ${last_lines_in_problem_reads_reformatted[@]}; do
    first_line_to_delete_in_seq=$((${line_num} - 3))
    second_line_to_delete_in_seq=$((${line_num} - 2))
    third_line_to_delete_in_seq=$((${line_num} - 1))
    fourth_line_to_delete_in_seq=${line_num}

    all_lines_to_delete[${index_for_all_lines_to_delete}]=${first_line_to_delete_in_seq}
    all_lines_to_delete[$((${index_for_all_lines_to_delete} + 1))]=${second_line_to_delete_in_seq}
    all_lines_to_delete[$((${index_for_all_lines_to_delete} + 2))]=${third_line_to_delete_in_seq}
    all_lines_to_delete[$((${index_for_all_lines_to_delete} + 3))]=${fourth_line_to_delete_in_seq}
    

    
    # sed "${starting_line_num_to_delete},${line_num}d" ${fastq_to_check} > polished.fq

    # Prepare for next iteration
    index_for_all_lines_to_delete=$((${index_for_all_lines_to_delete} + 4))
    # echo ${index_for_all_lines_to_delete}
done
# echo "${all_lines_to_delete[@]}d;"
# https://stackoverflow.com/a/48744059/8423001
# https://stackoverflow.com/a/26569006/8423001
# https://stackoverflow.com/a/15978536/8423001
# https://stackoverflow.com/a/26727351/8423001
# The first sed command replaces spaces and the first newline with "d;".
# It is used to format the 2nd sed command which does the actual deleting.
echo ${all_lines_to_delete[@]} | sed "s/\ /d;/g;s/$/d;/" \
    | xargs -I z sed z ${fastq_to_check} > ../../testing_data/sequencing_reads/polished.fq


