#!/bin/bash

fastqc \
    --outdir . \
    --extract \
    --threads 4 \
    --dir . \
    CGT3068_1.fq.gz