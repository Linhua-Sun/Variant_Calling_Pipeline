#!/bin/bash
#purpose: this script is used to call snp and indel applying two methods (GATK + samtools/bcftools) from DNA resequenceing DATA)
#author: Linhua Sun
#time: Fri Feb 26 08:20:07 CST 2016

#===========================================================================================================
#assign workspace                                                                                           
THREADS=10                                                                                                         
WORKSPACE=`pwd`                                                                                             
                                                                                                            
#assign reference genome of the species                                                                     
                                                                                                            
REF="${WORKSPACE}/REFERENCE_GENOME/IRGSP-1.0_genome.fasta"                                                  
                                                                                                            
#assign the location of softwares                                                                           
                                                                                                            
#GATK location                                                                                              
GATK="${WORKSPACE}/SOFTWARE/GATK/GenomeAnalysisTK.jar"                                                      
#PICARD location                                                                                            
PICARD_LOC="${WORKSPACE}/SOFTWARE/picard-tools-1.119"                                                       
                                                                                                            
                                                                                                            
                                                                                                            
#Assign input data variables (eg. ERR009626; `sh script.sh ERR009626`;$1)                             
                                                                                                            
SAMPLE="$WORKSPACE/RAW_DATA/$1"                                                                             
                                                                                                            
OUTPUT="$WORKSPACE/GATKDIR/$1"                                                                              
                                                                                                            
TEMP="$WORKSPACE/TEMP"


LOG=${OUTPUT}/${1}_LOG

if [ ! -d $LOG ]
		then mkdir -p $LOG
fi                                                                                    
                                                                                                            
if [ ! -d $OUTPUT ]                                                                                         
		then mkdir -p $OUTPUT                                                                               
fi                                                                                                          
                                                                                                            
if [ ! -d $TEMP ]                                                                                           
		then mkdir -p $WORKSPACE/TEMP                                                                       
fi                                                                                                          
#===========================================================================================================



#==============ONLY for the 1st time: Creating the different kinds index files==============================

#1.GATK create .dict file from $REF 
#do not delete &

#java -Xmx4G \
#    -jar $PICARD_LOC/CreateSequenceDictionary.jar \
#    R=$REF \
#    O=$WORKSPACE/REFERENCE_GENOME/IRGSP-1.0_genome.dict &

#2.samtools create .fai file from $REF

#samtools faidx $REF &

#3.bwa index; only for the first time

#bwa index $REF

#===========================================================================================================



#------------------------------------2).mapping-------------------------------------------------------------

SECONDS=0
#bwa -R change
bwa mem -M -t 10 -R "@RG\tID:"$1"\tSM:"$1"\tLB:"$1"\tPL:illumina\tPU:run" $REF $SAMPLE/${1}_1.fastq.gz $SAMPLE/${1}_2.fastq.gz > $OUTPUT/${1}.sam

echo "$SECONDS  bwa mem" > $OUTPUT/${1}_time.txt

SECONDS=0

#-------------------trans format

samtools fixmate -O bam $OUTPUT/${1}.sam $OUTPUT/${1}.bam

#-------------------sort

samtools sort -O bam -o $OUTPUT/${1}_sorted.bam -T $TEMP $OUTPUT/${1}.bam

#-------------------MARK duplicates

#Run the following Picard command to mark duplicates AND build index:

java -Xmx4G -jar $PICARD_LOC/MarkDuplicates.jar \
	INPUT=$OUTPUT/${1}_sorted.bam \
	OUTPUT=$OUTPUT/${1}_sorted_dedup_reads.bam \
	METRICS_FILE=$OUTPUT/${1}_metrics.txt \
	ASSUME_SORTED=TRUE

java -Xmx4G -jar $PICARD_LOC/BuildBamIndex.jar \
	INPUT=$OUTPUT/${1}_sorted_dedup_reads.bam

echo "$SECONDS  sort trans mark" >> $OUTPUT/${1}_time.txt

#----------------------------------3).Improvement-----------------------------------------------------------

#GATK Realigner. nct does not support here.
SECONDS=0
java -Xmx4G -jar $GATK \
	-nt $2 \
	-T RealignerTargetCreator \
	-R $REF \
	-I $OUTPUT/${1}_sorted_dedup_reads.bam \
	-o $OUTPUT/${1}_sorted_dedup_reads.bam.intervals

java -Xmx4G -jar $GATK \
	-T IndelRealigner \
	-R $REF \
	-I $OUTPUT/${1}_sorted_dedup_reads.bam \
	-targetIntervals $OUTPUT/${1}_sorted_dedup_reads.bam.intervals \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam

echo "$SECONDS  realign." >> $OUTPUT/${1}_time.txt

#samtools index outputted bam files

samtools index $OUTPUT/${1}_sorted_dedup_reads_realigned.bam

#----------------------------------3).known sites pre------------------------------------------------------
#1.HaplotypeCaller snp calling

#HaplotypeCaller ??? 
SECONDS=0
java -Xmx4G -jar $GATK \
	-R $REF \
	-T HaplotypeCaller \
	-nct $2 \
	--genotyping_mode DISCOVERY \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.gatk.raw1.vcf \
	--read_filter BadCigar \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10

echo "$SECONDS  HaplotypeCaller " >> $OUTPUT/${1}_time.txt

#2.samtools snp calling

SECONDS=0

samtools mpileup -ugf $REF $OUTPUT/${1}_sorted_dedup_reads_realigned.bam | bcftools call -vmO v -o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.samtools.raw1.vcf

echo "$SECONDS  samtools snp calling" >> $OUTPUT/${1}_time.txt

#3.merge vcf

#SelectVariants ???
SECONDS=0
java -Xmx20g -jar $GATK \
	-R $REF \
	-T SelectVariants \
	--variant $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.gatk.raw1.vcf \
	--concordance $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.samtools.raw1.vcf \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.vcf

MEANQUAL=`awk '{ if ($1 !~ /#/) { total += $6; count++ } } END { print total/count }' $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.vcf `

#VariantFiltration ???
java -Xmx4G -jar $GATK \
	-R $REF \
	-T VariantFiltration \
	--filterExpression "QD < 2.0 || ReadPosRankSum < -8.0 || FS > 60.0 || QUAL < $MEANQUAL" \
	--filterName LowQualFilter \
	--variant $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.vcf \
	--missingValuesInExpressionsShouldEvaluateAsFailing \
	--logging_level ERROR \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.flt1.vcf

#grep "LowQual" and "Filter"

grep -v "LowQual" \
$OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.flt1.vcf > \
$OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.flt1.confidence.rawlow.vcf

grep -v "Filter" \
$OUTPUT/${1}_sorted_dedup_reads_realigned.bam.concordance.raw.flt1.confidence.rawlow.vcf > \
$OUTPUT/${1}_final_raw.vcf

echo "$SECONDS  generate known vcf" >> $OUTPUT/${1}_time.txt

#===========================================================================================================

#Realigner
SECONDS=0
java -Xmx4G -jar $GATK \
	-nt $2 \
	-R $REF \
	-T RealignerTargetCreator \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam \
	-known $OUTPUT/${1}_final_raw.vcf \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.intervals

java -Xmx4G -jar $GATK \
	-R $REF \
	-T IndelRealigner \
	-targetIntervals $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.intervals \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam \
	-known $OUTPUT/${1}_final_raw.vcf \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam
echo "$SECONDS  Realigner_again " >> $OUTPUT/${1}_time.txt

#BaseRecalibrator
SECONDS=0

java -Xmx4G -jar $GATK \
	-nct $2 \
	-R $REF \
	-T BaseRecalibrator \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam.recal_data.grp \
	-knownSites $OUTPUT/${1}_final_raw.vcf

java -Xmx4G -jar $GATK \
	-nct $2 \
	-T BaseRecalibrator \
	-R $REF \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam \
	-BQSR $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam.recal_data.grp \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam.recal_data.grp2 \
	-knownSites $OUTPUT/${1}_final_raw.vcf

#AnalyzeCovariates
java -jar $GATK \
	-T AnalyzeCovariates \
	-R $REF \
	-before $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam.recal_data.grp \
	-after $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam.recal_data.grp2 \
	-csv $OUTPUT/${1}_recal.grp.csv \
	-plots $OUTPUT/${1}_recal.grp.pdf

#PrintReads trans it into bam file
java -Xmx4G -jar $GATK \
	-nct $2 \
	-R $REF \
	-T PrintReads \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam \
	-o $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.recal.bam \
	-BQSR $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.realn.bam.recal_data.grp

echo "$SECONDS  BaseRecalibrator " >> $OUTPUT/${1}_time.txt

#HaplotypeCaller again:

SECONDS=0

java -Xmx4G -jar $GATK \
	-R $REF \
	-nct $2 \
	-T HaplotypeCaller \
	--genotyping_mode DISCOVERY \
	-I $OUTPUT/${1}_sorted_dedup_reads_realigned.bam.recal.bam \
	-o $OUTPUT/${1}_gatk.raw2.vcf \
	--read_filter BadCigar \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10

echo "$SECONDS  HaplotypeCaller again " >> $OUTPUT/${1}_time.txt

#SelectVariants
SECONDS=0

java -Xmx4G -jar $GATK \
	-R $REF \
	-T SelectVariants \
	-V $OUTPUT/${1}_gatk.raw2.vcf \
	-selectType SNP \
	-o $OUTPUT/${1}_gatk.raw_snp.vcf

java -Xmx4G -jar $GATK \
	-R $REF \
	-T SelectVariants \
	-V $OUTPUT/${1}_gatk.raw2.vcf \
	-selectType INDEL \
	-o $OUTPUT/${1}_gatk.raw_indel.vcf

#VariantFiltration SNP

MEANQUAL1=$(awk '{ if ($1 !~ /#/) { total += $6; count++ } } END { print total/count }' $OUTPUT/${1}_gatk.raw_snp.vcf)

java -Xmx4G -jar $GATK \
	-R $REF \
	-T VariantFiltration \
	--filterExpression "MQ < 40.0 || MQRankSum < -12.5 || QD < 2.0 || ReadPosRankSum < -8.0 || FS > 60.0 || QUAL < $MEANQUAL1" \
	--filterName LowQualFilter \
	--variant $OUTPUT/${1}_gatk.raw_snp.vcf \
	--missingValuesInExpressionsShouldEvaluateAsFailing \
	--logging_level ERROR \
	-o $OUTPUT/${1}.concordance.flt_snp.vcf

grep -v "Filter" \
	$OUTPUT/${1}.concordance.flt_snp.vcf > \
	$OUTPUT/${1}.concordance.flt_snp.final.vcf

grep -v "LowQual" \
	$OUTPUT/${1}.concordance.flt_snp.final.vcf > \
	$OUTPUT/${1}.concordance.flt_snp.fina2.vcf

#VariantFiltration INDEL

MEANQUAL2=`awk '{ if ($1 !~ /#/) { total += $6; count++ } } END { print total/count }' $OUTPUT/${1}_gatk.raw_indel.vcf`

java -Xmx4G -jar $GATK \
	-R $REF \
	-T VariantFiltration \
	--filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 ||  FS > 200.0 || QUAL < $MEANQUAL2" \
	--filterName LowQualFilter \
	--variant $OUTPUT/${1}_gatk.raw_indel.vcf \
	--missingValuesInExpressionsShouldEvaluateAsFailing \
	--logging_level ERROR \
	-o $OUTPUT/${1}.concordance.flt_indel.vcf

grep -v "Filter" \
	$OUTPUT/${1}.concordance.flt_indel.vcf > \
	$OUTPUT/${1}.concordance.flt_indel.final.vcf

grep -v "LowQual" \
	$OUTPUT/${1}.concordance.flt_indel.final.vcf > \
	$OUTPUT/${1}.concordance.flt_indel.fina2.vcf

echo "$SECONDS  finished filter " >> $OUTPUT/${1}_time.txt