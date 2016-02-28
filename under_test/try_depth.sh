#bin/bash
#java -jar /sdd2/users/linhua/11.1.2016/SOFTWARE/GATK/GenomeAnalysisTK.jar \
#	-T DepthOfCoverage \
#	-R /sdd2/users/linhua/11.1.2016/REFERENCE_GENOME/IRGSP-1.0_genome.fasta \
#	-o file_name_base \
#	-I GATKDIR/ERR037209/ERR037209_sorted_dedup_reads_realigned.bam
#	#[-geneList refSeq.sorted.txt] \
#	#[-pt readgroup] \
#	#[-ct 4 -ct 6 -ct 10] \
#	#[-L my_capture_genes.interval_list]	
#	#/sdd2/users/linhua/11.1.2016/SOFTWARE/GATK/GenomeAnalysisTK.jar
#	#/sdd2/users/linhua/11.1.2016/REFERENCE_GENOME/IRGSP-1.0_genome.fasta


java -Xmx4g -jar /sdd2/users/linhua/11.1.2016/SOFTWARE/GATK/GenomeAnalysisTK.jar \
	-R /sdd2/users/linhua/11.1.2016/REFERENCE_GENOME/IRGSP-1.0_genome.fasta \
	-T DepthOfCoverage \
	--omitDepthOutputAtEachBase \
	--omitIntervalStatistics \
	--omitLocusTable \
	-ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 \
	--nBins 99 --start 1 --stop 100 \
	-I GATKDIR/ERR037209/ERR037209_sorted_dedup_reads_realigned.bam \
	-o sample.depth
