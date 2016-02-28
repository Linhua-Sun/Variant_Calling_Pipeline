#!/bin/bash
#===========================================================================================================
#assign variables

WORKSPACE=$(pwd)

#assign reference genome of the species

REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#GATK location
GATK="${WORKSPACE}/SOFTWARE/GATK/GenomeAnalysisTK.jar"

OUTPUT2="${WORKSPACE}/GATKDIR/GATK-GVCF-test"

if [ ! -d $OUTPUT2 ]
		then mkdir -p $OUTPUT2
fi

# Perform joint genotyping on gVCF files produced by HaplotypeCaller
#find and xargs echos the standard expression GATK need. http://www.51testing.com/html/32/15077732-1554471.html

GVCFS=`find ./ -name "*_output_raw_snps_indels.g.vcf" -print | \
    xargs -I GVCF_FILE echo -n "-V GVCF_FILE "`

java -Xmx10G -jar $GATK \
	-R $REF \
	-T GenotypeGVCFs \
	-nt 5 \
	-o ${OUTPUT2}/GATK_gt_raw.vcf \
	${GVCFS} \
	> ${OUTPUT2}/GATK_GVCF.log 2>&1
