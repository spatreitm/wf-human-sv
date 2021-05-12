#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage(){
    log.info """
wf-human-sv

Usage:
    nextflow run epi2melabs/wf-human-sv [options]

Script Options:
    --fastq             DIR     FASTQ file (recommended)
    --reference         FILE    Reference FASTA file (required)
    --out_dir           DIR     Path for output (default: $params.out_dir)
    --bam               DIR     BAM file (as alternative to FASTQ)
    --bed               FILE    Target bed file
    --eval              BOOL    Perform benchmarking against truthset
    --help
    
"""
}


