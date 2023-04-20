#!/bin/bash
#Joel Markus
#Run amrfinder

#activate the environment
conda activate amrfinder

#Run using input fasta, output to STDOUT
amrfinder --protein ${1} --threads 10 --plus > ${2} 

#deactivate the environment
conda deactivate
