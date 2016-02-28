#!/bin/bash

#Linhua Sun

#Fri Jan 29 10:56:09 CST 2016

#Assign workspace

WORKSPACE="./"

#assign reference genome of rice

REF="./REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares(speedseq):

speedseq="/sdd2/users/linhua/11.1.2016/SOFTWARE/speedseq/bin/speedseq"

#assign the location of files:

SAMPLE="$WORKSPACE/RAW_DATA/$1"

OUTPUT="$WORKSPACE/SPEEDDIR/$1"

TEMP="$WORKSPACE/TEMP"

if [ ! -d $OUTPUT ]
		then mkdir -p $OUTPUT
fi

if [ ! -d $TEMP ]
		then mkdir -p $TEMP
fi

#test
#	RAW_DATA/SMALL_DATA_1:
#	SMALL_DATA_1_1.fastq.gz  SMALL_DATA_1_2.fastq.gz


#Call variants on multiple samples
#Use speedseq align to produce sorted, duplicate-marked, BAM alignments for each sample.

#$speedseq align \
#	-t 8 \
#	-o $OUTPUT/$1 \
#	-R "@RG\tID:$1\tSM:$1\tLB:$1\tPL:illumina\tPU:run" \
#	$REF \
#	$SAMPLE/${1}_1.fastq.gz \
#	$SAMPLE/${1}_2.fastq.gz 

#Use speedseq var to call SNVs and indels on multiple samples.

$speedseq var \
	-o $OUTPUT/$1 \
	-t 3 \
	$REF \
	$OUTPUT/$1.bam
#
#-w annotations/ceph18.b37.include.2014-01-15.bed \