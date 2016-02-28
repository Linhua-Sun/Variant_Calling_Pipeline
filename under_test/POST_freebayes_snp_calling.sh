#!/bin/bash

#freebayes snp calling

#===========================================================================================================
#assign variables

WORKSPACE=$(pwd)

#assign reference genome of the species

REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

OUTPUT="${WORKSPACE}/freebayesDIR/"

if [ ! -d $OUTPUT ]
		then mkdir -p $OUTPUT
fi

freebayes -L ready_bams_path.txt --fasta-reference $REF > ready_bams_freebayes.vcf

#-L --bam-list FILE :A file containing a list of BAM files to be analyzed.
