#!/bin/bash
#Assign workspace

WORKSPACE=`pwd`

#assign reference genome of rice

REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares(speedseq):

speedseq="/sdd2/users/linhua/11.1.2016/SOFTWARE/speedseq/bin/speedseq"

$speedseq var \
	-o speedseq_fb_test \
	-t 8 \
	-T TEMP \
	$REF \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037208/ERR037208_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037207/ERR037207_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037206/ERR037206_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037205/ERR037205_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037204/ERR037204_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037203/ERR037203_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037202/ERR037202_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037201/ERR037201_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037200/ERR037200_sorted_dedup_reads_realigned.bam.recal.bam \
	/sdd2/users/linhua/11.1.2016/GATKDIR/ERR037209/ERR037209_sorted_dedup_reads_realigned.bam.recal.bam

	