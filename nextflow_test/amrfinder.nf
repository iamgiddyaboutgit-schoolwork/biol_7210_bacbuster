#!/usr/bin/env nextflow

// Define the process to run the amrfinder script
process run_amrfinder {
    input:
    file input_fasta

    output:
    file "${params.output}"

    script:
    """
    #!/bin/bash
    #Joel Markus
    #Run amrfinder

    #activate the environment
    conda activate amrfinderplus

    #Run using input fasta, output to STDOUT
    amrfinder --protein ${input_fasta} --threads ${params.threads} --plus > ${params.output}

    #deactivate the environment
    conda deactivate
    """
}

// Define the input parameters
params.input = "./fasta_n/CGT3027_S_genes.fa"
params.output = "amrfindertest_output.txt"
params.threads = 10

// Define the workflow
workflow {
    // Define the input fasta file
    input_fasta = file(params.input)

    // Run the amrfinder process
    run_amrfinder(input_fasta)
        .params(params)
        .withName('amrfinder')
        .withReport(true)
}
