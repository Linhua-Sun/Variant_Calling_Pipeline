#!/bin/bash
#purpose: this script is used to call snp and indel applying two methods (GATK + samtools/bcftools) from DNA resequenceing DATA)
#author: Linhua Sun
#time: Sun Jan 31 21:46:30 CST 2016

#===========================================================================================================
#assign workspace

WORKSPACE="./"

#assign reference genome of the species

REF="./REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares
#GATK location
GATK="./SOFTWARE/GATK/GenomeAnalysisTK.jar"
#PICARD location
PICARD_LOC="./SOFTWARE/picard-tools-1.119"

#samtools version (be careful about it!)
#Program: samtools (Tools for alignments in the SAM format)
#Version: 1.3 (using htslib 1.3)
#Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
#Version: 1.3 (using htslib 1.3)

#Assign input data variables (eg. ERR009626; `sh script.sh ERR009626 10`;$1;$2)

SAMPLE="$WORKSPACE/RAW_DATA/$1"

OUTPUT="$WORKSPACE/GATKDIR/$1"

TEMP="$WORKSPACE/TEMP"

if [ ! -d $OUTPUT ]
		then mkdir -p $OUTPUT
fi

if [ ! -d $TEMP ]
		then mkdir -p $WORKSPACE/TEMP
fi
bwa mem -M -t 2 -R "@RG\tID:"$1"\tSM:"$1"\tLB:"$1"\tPL:illumina\tPU:run" $REF $SAMPLE/${1}_1.fastq.gz $SAMPLE/${1}_2.fastq.gz > $OUTPUT/${1}.sam