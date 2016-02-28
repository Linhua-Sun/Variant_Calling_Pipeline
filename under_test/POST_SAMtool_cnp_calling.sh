#!/bin/bash
#===========================================================================================================
#assign variables

WORKSPACE=$(pwd)

#assign reference genome of the species

REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

OUTPUT="${WORKSPACE}/SAMtoolsDIR/"

if [ ! -d $OUTPUT ]
		then mkdir -p $OUTPUT
fi
samtools mpileup -ugf $REF -b ready_bams_path.txt| bcftools call -vmO z -o $OUTPUT/ready_bams_SAM.vcf.gz

# -b, --bam-list FILE  List of input BAM files, one file per line [null]
