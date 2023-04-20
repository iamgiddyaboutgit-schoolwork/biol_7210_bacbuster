#!/usr/bin/env bash

# Cleans up nextflow logs and work directory.

rm -r -f -v ".nextflow"
rm -r -f -v "work"
rm -f -v ".nextflow.log"*
