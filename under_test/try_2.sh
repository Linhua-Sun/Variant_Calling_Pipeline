#!/bin/bash
#Assign workspace

WORKSPACE="./"

#assign reference genome of rice

REF="./REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares(speedseq):

speedseq="/sdd2/users/linhua/11.1.2016/SOFTWARE/speedseq/bin/speedseq"

#assign the location of files:

OUTPUT="$WORKSPACE/SPEEDDIR/"

if [ ! -d $OUTPUT ]
		then mkdir -p $OUTPUT
fi

$speedseq var \
	-o $OUTPUT \
	-t 5 \
	$REF \
	./SPEEDDIR/ERR037208/ERR037208.bam \
	./SPEEDDIR/ERR037207/ERR037207.bam \
	./SPEEDDIR/ERR037206/ERR037206.bam \
	./SPEEDDIR/ERR037205/ERR037205.bam \
	./SPEEDDIR/ERR037204/ERR037204.bam \
	./SPEEDDIR/ERR037203/ERR037203.bam \
	./SPEEDDIR/ERR037202/ERR037202.bam \
	./SPEEDDIR/ERR037201/ERR037201.bam \
	./SPEEDDIR/ERR037200/ERR037200.bam \
	./SPEEDDIR/ERR037209/ERR037209.bam

#-w forbided

