#!/bin/bash

#test the function of GenotypeConcordance
#https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeConcordance.php
#===========================================================================================================
#assign variables

WORKSPACE=/sdd2/users/linhua/11.1.2016

#assign reference genome of the species

REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares

#GATK location
GATK="${WORKSPACE}/SOFTWARE/GATK/GenomeAnalysisTK.jar"

java -Xmx10G -jar $GATK \
	-T GenotypeConcordance \
	-R $REF \
	-eval /sdd2/users/linhua/11.1.2016/vcf_gz_stats/GATK_gt_raw.vcf.gz \
	-comp /sdd2/users/linhua/11.1.2016/vcf_gz_stats/ready_bams_SAM.vcf.gz \
	-o output.grp
