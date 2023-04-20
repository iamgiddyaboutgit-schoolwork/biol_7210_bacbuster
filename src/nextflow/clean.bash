#!/usr/bin/env bash

# Cleans up nextflow logs and work directory.

rm -r -f -v ".nextflow"
mv --target-directory="../" work/conda
rm -r -f -v "work"
rm -f -v ".nextflow.log"*
mkdir -p work
mv --target-directory="work" ../conda
