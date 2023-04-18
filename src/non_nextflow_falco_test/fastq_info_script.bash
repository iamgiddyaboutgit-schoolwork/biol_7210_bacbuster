#!/usr/bin/env bash
fastq_to_check="../../testing_data/sequencing_reads/test_fastq_3.fq"
num_lines_in_file=$(wc -l ${fastq_to_check} | cut -f 1 -d " ")
num_reads_in_file=$((${num_lines_in_file} / 4))

#####################################################################
# Check for specific file format problems.
# Note that fastq_info will terminate upon finding the first error.
# Using the trick here,
# https://stackoverflow.com/a/2342841/8423001
# get the line with the problem.
 
# last_lines_in_problem_reads gives the last line in the 4-line sequence with the problem.
# Problems that we use a regex to find are "duplicated sequence" and "sequence and 
# quality score lengths don't match".
# Setup before loop:
problem_counter=0

last_lines_in_problem_reads[${problem_counter}]=$(fastq_info ${fastq_to_check} 2>&1 > /dev/null \
    | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)

if  [[ -n "${last_lines_in_problem_reads[0]}" ]]
then
    # last_lines_in_problem_reads[0] is of non-zero length, i.e.
    # we found a problem.
    line_after_problem=$((${last_lines_in_problem_reads[${problem_counter}]} + 1))
    ((problem_counter+=1))
else
    # Stop script; there were no problems.
    exit 0
fi

echo ${line_after_problem}
# Iteratively check for additional problems if there was a 
# previous problem.
while (( "${line_after_problem}" < "${num_lines_in_file}" ))
do
    echo "line_after_problem: ${line_after_problem}"
    echo "num_lines_in_file: ${num_lines_in_file}"
    echo $(( "${line_after_problem}" < "${num_lines_in_file}" ))
    last_lines_in_problem_reads[${problem_counter}]=$(tail -n +${line_after_problem} ${fastq_to_check} \
        | fastq_info - 2>&1 > /dev/null \
        | grep -P --max-count=1 --only-matching "(?<=line\s)[0-9]+(?=:\s((duplicated\sseq)|(sequence\sand\squality)))" -)
    
    echo "array: ${last_lines_in_problem_reads[@]}"
    if  [[ -n "${last_lines_in_problem_reads[${problem_counter}]}" ]]
    then 
        echo "We found another problem."
        # We found another problem.
        line_after_problem=$((${last_lines_in_problem_reads[${problem_counter}]} + 1))
        ((problem_counter+=1))
    else
        # Stop script; there were no further problems (of the kind we are looking for)
        exit 0
    fi
    
done

echo ${last_lines_in_problem_reads}
# # We assume that we want to delete all 4 lines with the problem in the .fastq file.

# # A loop can be used to get a running sum of the elements of 
# # last_lines_in_problem_reads.  From this, we can infer
# # the line numbers in the original file.
# # The idea is that if we want to do all of our deleting of lines
# # in the file in memory, then we need to know the actual line
# # numbers that need to be deleted ahead of time.  
# # We can easily get an array of line numbers using fastq_info
# # that are based on lines since the last recursive call.
# # To change this array to absolute line numbers, we just
# # need to perform some simple arithmetic.
# # Declare an array.
# declare -a last_lines_in_problem_reads_reformatted 

# num_last_lines_in_problem_reads=${#last_lines_in_problem_reads[@]}

# if  [[ -n "${last_lines_in_problem_reads[0]}" ]]; then
#     # last_lines_in_problem_reads[0] is of non-zero length
#     last_lines_in_problem_reads_reformatted[0]=${last_lines_in_problem_reads[0]}
# fi

# j=1
# # After this loop, last_lines_in_problem_reads_reformatted will hold
# # line numbers positioned within the original file.
# for (( i=1 ; i<${num_last_lines_in_problem_reads} ; i++ )); do
#     last_lines_in_problem_reads_reformatted[${j}]=$((${last_lines_in_problem_reads_reformatted[$((${j} - 1))]} + ${last_lines_in_problem_reads[i]}))
#     j=$((${j}+1))
    
# done

# declare -a all_lines_to_delete
# index_for_all_lines_to_delete=0
# for line_num in ${last_lines_in_problem_reads_reformatted[@]}; do
#     first_line_to_delete_in_seq=$((${line_num} - 3))
#     second_line_to_delete_in_seq=$((${line_num} - 2))
#     third_line_to_delete_in_seq=$((${line_num} - 1))
#     fourth_line_to_delete_in_seq=${line_num}

#     all_lines_to_delete[${index_for_all_lines_to_delete}]=${first_line_to_delete_in_seq}
#     all_lines_to_delete[$((${index_for_all_lines_to_delete} + 1))]=${second_line_to_delete_in_seq}
#     all_lines_to_delete[$((${index_for_all_lines_to_delete} + 2))]=${third_line_to_delete_in_seq}
#     all_lines_to_delete[$((${index_for_all_lines_to_delete} + 3))]=${fourth_line_to_delete_in_seq}

#     # Prepare for next iteration
#     index_for_all_lines_to_delete=$((${index_for_all_lines_to_delete} + 4))
# done

# # https://stackoverflow.com/a/48744059/8423001
# # https://stackoverflow.com/a/26569006/8423001
# # https://stackoverflow.com/a/15978536/8423001
# # https://stackoverflow.com/a/26727351/8423001
# # The first sed command replaces spaces and the first newline with "d;".
# # It is used to format the 2nd sed command which does the actual deleting.
# echo ${all_lines_to_delete[@]} | sed "s/\ /d;/g;s/$/d;/" \
#     | xargs -I z sed z ${fastq_to_check} > ../../testing_data/sequencing_reads/polished.fq


