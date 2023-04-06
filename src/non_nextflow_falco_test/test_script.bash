#!/bin/bash

falco \
	--outdir . \
        --threads 1 \
        -subsample 1000 \
        -skip-data \
        -skip-report \
        CGT3068_1.fq.gz
