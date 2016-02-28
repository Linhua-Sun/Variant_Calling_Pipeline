#!/bin/bash
#===========================================================================================================
#assign variables

WORKSPACE=$(pwd)

#assign reference genome of the species

REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"

#assign the location of softwares

#GATK location
GATK="${WORKSPACE}/SOFTWARE/GATK/GenomeAnalysisTK.jar"

#PICARD location
PICARD_LOC="${WORKSPACE}/SOFTWARE/picard-tools-1.119"

#samtools version (be careful about it!)
#Program: samtools (Tools for alignments in the SAM format)
#Version: 1.3 (using htslib 1.3)
#Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
#Version: 1.3 (using htslib 1.3)

#Assign input data variables (eg. ERR009626; `sh script.sh ERR009626 10`;$1;$2)

SAMPLE="${WORKSPACE}/RAW_DATA/$1"

OUTPUT="${WORKSPACE}/GATKDIR/$1"

TEMP="${WORKSPACE}/TEMP"

if [ ! -d $OUTPUT ]
		then mkdir -p $OUTPUT
fi

if [ ! -d $TEMP ]
		then mkdir -p $WORKSPACE/TEMP
fi

#input bam file: ${OUTPUT}/${1}_sorted_dedup_reads_realigned.bam.recal.bam
#onput gvcf file: ${OUTPUT}/${1}_output_raw_snps_indels.g.vcf

java -jar $GATK \
	-nct $2 \
	-R $REF \
	-T HaplotypeCaller \
	-I ${OUTPUT}/${1}_sorted_dedup_reads_realigned.bam.recal.bam \
	--emitRefConfidence GVCF \
	-o ${OUTPUT}/${1}_output_raw_snps_indels.g.vcf

#1. Variant calling
#Run the HaplotypeCaller on each sample's BAM file(s) (if a sample's data is spread over more than one BAM, then pass them all in together) to create single-sample gVCFs, with the option -emitRefConfidence GVCF, and using the .g.vcf extension for the output file.
#
#
#2. Optional data aggregation step
#If you have more than a few hundred samples, run CombineGVCFs on batches of ~200 gVCFs to hierarchically merge them into a single gVCF. This will make the next step more tractable.
#
#3. Joint genotyping
#Take the outputs from step 2 (or step 1 if dealing with fewer samples) and run GenotypeGVCFs on all of them together to create the raw SNP and indel VCFs that are usually emitted by the callers.
#
#Genotype across samples:
#java -Xmx64g -jar ~/Programs/GenomeAnalysisTK.jar -R /reference.fasta -T GenotypeGVCFs -o woth300-rawsnps.vcf -nt 24 --variant sample1.g.vcf --variant sample2.g.vcf ... ... --variant sample 296.g.vcf




