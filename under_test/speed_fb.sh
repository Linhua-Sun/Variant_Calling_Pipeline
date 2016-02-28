#!/bin/bash
#Assign workspace

WORKSPACE="./"

#assign reference genome of rice

REF="./REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares(speedseq):

speedseq="/sdd2/users/linhua/11.1.2016/SOFTWARE/speedseq/bin/speedseq"

$speedseq var \
	-o FB-test \
	-t 10 \
	$REF \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037208/ERR037208.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037207/ERR037207.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037206/ERR037206.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037205/ERR037205.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037204/ERR037204.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037203/ERR037203.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037202/ERR037202.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037201/ERR037201.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037200/ERR037200.bam \
	/sdd2/users/linhua/11.1.2016/SPEEDDIR/ERR037209/ERR037209.bam

#-w forbided