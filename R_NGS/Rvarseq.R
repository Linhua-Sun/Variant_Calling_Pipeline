### R code from vignette source 'Rvarseq.Rnw'

###################################################
### code chunk number 1: Rvarseq.Rnw:265-266 (eval = FALSE)
###################################################
## download.file("http://biocluster.ucr.edu/~nkatiyar/Rvarseq_workshop_2014/Rvarseq.zip", "Rvarseq.zip")


###################################################
### code chunk number 2: Rvarseq.Rnw:268-270
###################################################
targets <- read.delim("./data/targets.txt")
targets


###################################################
### code chunk number 3: Rvarseq.Rnw:284-287 (eval = FALSE)
###################################################
## library(modules); library(Rsamtools)
## moduleload("bwa/0.7.10") # loads BWA version 0.7.10 from module system
## system("bwa index -a bwtsw ./data/tair10chr.fasta") # Indexes reference genome; required for GATK


###################################################
### code chunk number 4: Rvarseq.Rnw:291-297 (eval = FALSE)
###################################################
## dir.create("results") # Note: all output data will be written to results directory
## moduleload("bwa/0.7.10")
## for(i in seq(along=targets[,1])) {
##   system(paste("bwa mem -M -R '@RG\\tID:group1\\tSM:sample1\\tPL:illumina\\tLB:lib1\\tPU:unit1'", " ./data/tair10chr.fasta", " ./data/",targets$FileName[i], " >", " ./results/",targets$FileName[i], ".sam", sep=""))
##     asBam(file=paste("./results/", targets$FileName[i], ".sam", sep=""), destination=paste("./results/", targets$FileName[i], sep=""), overwrite=TRUE, indexDestination=TRUE)
##     }


###################################################
### code chunk number 5: Rvarseq.Rnw:308-311 (eval = FALSE)
###################################################
## library(gmapR); library(rtracklayer)
## fastaFile <- FastaFile(paste(getwd(), "/data/tair10chr.fasta", sep="")) # Needs to be full path!
## gmapGenome <- GmapGenome(fastaFile, directory="data", name="gmap_tair10chr/", create=TRUE)


###################################################
### code chunk number 6: Rvarseq.Rnw:317-324 (eval = FALSE)
###################################################
## gmapGenome <- GmapGenome(fastaFile, directory="data", name="gmap_tair10chr/", create=FALSE)
## # To regenerate gmapGenome object, set 'create=FALSE'.
## param <- GsnapParam(genome=gmapGenome, unique_only = TRUE, molecule = "DNA", max_mismatches = 3)
## for(i in seq(along=targets[,1])) {
## output <- gsnap(input_a=paste("./data/", targets[i,1], sep=""), input_b=NULL, param,
## output=paste("results/gsnap_bam/", targets[i,1], sep=""))
## }


###################################################
### code chunk number 7: Rvarseq.Rnw:338-349
###################################################
library(VariantTools); library(gmapR)
gmapGenome <- GmapGenome(genome="gmap_tair10chr", directory="data")
tally.param <- TallyVariantsParam(gmapGenome, high_base_quality = 23L, indels = TRUE)
bfl <- BamFileList(paste("./results/", as.character(targets[,1]), ".bam", sep=""), index=character())
var <- callVariants(bfl[[1]], tally.param)
length(var)
var <- var[totalDepth(var) == altDepth(var) & totalDepth(var)>=5 & values(var)$n.read.pos >= 5] # Some arbitrary filter
length(var)
sampleNames(var) <- "bwa"
vcf <- asVCF(var)
writeVcf(vcf, "./results/varianttools.vcf", index = TRUE)


###################################################
### code chunk number 8: Rvarseq.Rnw:355-361 (eval = FALSE)
###################################################
## bfl <- BamFileList(paste("./results/gsnap_bam/", as.character(targets[,1]), ".sam", ".bam", sep=""), index=character())
## var_gsnap <- callVariants(bfl[[1]], tally.param)
## var_gsnap <- var_gsnap[totalDepth(var_gsnap) == altDepth(var_gsnap) & totalDepth(var_gsnap)>=5 & values(var_gsnap)$n.read.pos >= 5]
## sampleNames(var_gsnap) <- "gsnap"
## vcf_gsnap <- asVCF(var_gsnap)
## writeVcf(vcf_gsnap, "./results/varianttools_gnsap.vcf", index=TRUE)


###################################################
### code chunk number 9: Rvarseq.Rnw:371-372
###################################################
raw.variants <- tallyVariants(bfl[[1]], tally.param)


###################################################
### code chunk number 10: Rvarseq.Rnw:377-379
###################################################
qa.variants <- qaVariants(raw.variants)
softFilterMatrix(qa.variants)[1:2,]


###################################################
### code chunk number 11: Rvarseq.Rnw:384-386
###################################################
called.variants <- callVariants(qa.variants)
length(called.variants)


###################################################
### code chunk number 12: Rvarseq.Rnw:400-404
###################################################
library(VariantAnnotation)
vcf_imported <- readVcf("results/varianttools.vcf.bgz", "ATH1")
VRangesFromVCF <- as(vcf_imported, "VRanges")
VRangesFromVCF[1:4,]


###################################################
### code chunk number 13: Rvarseq.Rnw:427-438 (eval = FALSE)
###################################################
## library(modules)
## moduleload("samtools")
## dedup <- paste("samtools rmdup -S ", path(bfl[[1]]), " ", path(bfl[[1]]), "dedup.bam", sep="")
## system(dedup) # Removes PCR duplicates with identical read mappings!
## indexBam(file=paste(path(bfl[[1]]), "dedup.bam", sep=""))
## vcf1 <- paste("samtools mpileup -uf ./data/tair10chr.fasta ", path(bfl[[1]]), "dedup.bam",
##     " | bcftools view -bvcg  -> ./results/sambcf.raw.bcf", sep="")
## vcf2 <- paste("bcftools view ./results/sambcf.raw.bcf",
##     "| vcfutils.pl varFilter -D100 > ./results/sambcf.vcf")
## system(vcf1)
## system(vcf2)


###################################################
### code chunk number 14: Rvarseq.Rnw:447-458 (eval = FALSE)
###################################################
## library(modules)
## moduleload("java")
## system("java -jar /opt/picard/1.81/CreateSequenceDictionary.jar R=data/tair10chr.fasta O=data/tair10chr.dict")
## dir.create("results/gatktmp", recursive = TRUE)
## file.copy("gatk_runs.sh", "results/gatktmp/gatk_runs.sh")
## file.copy("results/SRR064154.fastq.bam", "results/gatktmp/myfile.fastq.bam")
## setwd("results/gatktmp")
## system("./gatk_runs.sh")
## file.copy("vargatk.recalibrated.filtered.vcf", "../gatk.vcf")
## setwd("../../")
## unlink("results/gatktmp/", recursive=TRUE, force=TRUE)


###################################################
### code chunk number 15: Rvarseq.Rnw:473-483
###################################################
library(VariantAnnotation)
vcfsam <- readVcf("results/sambcf.vcf", "ATH1")
vcfvt <- readVcf("results/varianttools.vcf.bgz", "ATH1")
vcfvt_gsnap <- readVcf("results/varianttools_gnsap.vcf.bgz", "ATH1")
vcfgatk <- readVcf("results/gatk.vcf", "ATH1")
vcfgatk <- vcfgatk[values(rowData(vcfgatk))$FILTER == "PASS"] # Uses GATK filters
methods <- list(BCF_BWA=names(rowData(vcfsam)), VariantTools_BWA=names(rowData(vcfvt)), VariantTools_GSNAP=names(rowData(vcfvt_gsnap)), GATK_BWA=names(rowData(vcfgatk)))
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
OLlist <- overLapper(setlist=methods, sep="_", type="vennsets")
counts <- sapply(OLlist$Venn_List, length); vennPlot(counts=counts, mymain="Variant Calling Methods")


###################################################
### code chunk number 16: Exercise 1 (eval = FALSE)
###################################################
## bfl <- BamFileList(paste("./results/", as.character(targets[,1]), ".bam", sep=""), index=character())
## vrl <- sapply(basename(names(bfl)), function(x) NULL, simplify=FALSE)
## for(i in seq(along=bfl)) {
##     var <- callVariants(bfl[[i]], tally.param)
##     var <- var[totalDepth(var) == altDepth(var) & totalDepth(var)>=5 & values(var)$n.read.pos >= 5]
##     vrl[[i]] <- paste(as.character(seqnames(var)), ":", start(var), "_", ref(var), "/", alt(var), sep="")
##     print(i)
##     }   
## source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
## OLlist <- overLapper(setlist=vrl, sep="_", type="vennsets")
## counts <- sapply(OLlist$Venn_List, length); vennPlot(counts=counts, mymain="Variant Calls for Four Samples")
## OLlist[[4]][15]


###################################################
### code chunk number 17: Rvarseq.Rnw:529-538
###################################################
library(GenomicFeatures)
chrominfo <- data.frame(chrom=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"), length=rep(10^5, 7), is_circular=rep(FALSE, 7))
txdb <- makeTranscriptDbFromGFF(file="data/TAIR10_GFF3_trunc.gff",
        format="gff3",
        dataSource="TAIR",
        chrominfo=chrominfo,
        species="Arabidopsis thaliana")
saveDb(txdb, file="./data/TAIR10.sqlite")
txdb <- loadDb("./data/TAIR10.sqlite")


###################################################
### code chunk number 18: Rvarseq.Rnw:541-544
###################################################
library(VariantAnnotation)
vcf <- readVcf("results/varianttools_gnsap.vcf.bgz", "ATH1")
seqlengths(vcf) <- seqlengths(txdb)[names(seqlengths(vcf))]; isCircular(vcf) <- isCircular(txdb)[names(seqlengths(vcf))]


###################################################
### code chunk number 19: Rvarseq.Rnw:547-549
###################################################
library(Rsamtools)
fa <- FaFile("data/tair10chr.fasta")


###################################################
### code chunk number 20: Rvarseq.Rnw:565-567
###################################################
vcf <- readVcf("results/sambcf.vcf", "ATH1")
seqlengths(vcf) <- seqlengths(txdb)[names(seqlengths(vcf))]; isCircular(vcf) <- isCircular(txdb)[names(seqlengths(vcf))]


###################################################
### code chunk number 21: Rvarseq.Rnw:577-579
###################################################
seqinfo(vcf)
genome(vcf)


###################################################
### code chunk number 22: Rvarseq.Rnw:590-594
###################################################
header(vcf)
meta(header(vcf))
info(header(vcf))[1:3,]
geno(header(vcf))[1:3,]


###################################################
### code chunk number 23: Rvarseq.Rnw:606-607
###################################################
rowData(vcf)[1:3,] 


###################################################
### code chunk number 24: Rvarseq.Rnw:610-611
###################################################
info(vcf)[1:3,1:6]


###################################################
### code chunk number 25: Rvarseq.Rnw:614-615
###################################################
alt(vcf)[1:3,]


###################################################
### code chunk number 26: Rvarseq.Rnw:628-633
###################################################
library(GenomicFeatures)
vcf <- readVcf(file="results/varianttools_gnsap.vcf.bgz", genome="ATH1")
seqlengths(vcf) <- seqlengths(txdb)[names(seqlengths(vcf))]; isCircular(vcf) <- isCircular(txdb)[names(seqlengths(vcf))]
rd <- rowData(vcf)
codvar <- locateVariants(rd, txdb, CodingVariants())


###################################################
### code chunk number 27: Rvarseq.Rnw:667-669
###################################################
allvar <- locateVariants(rd, txdb, AllVariants())
allvar[1:4]


###################################################
### code chunk number 28: Rvarseq.Rnw:674-677
###################################################
source("Rvarseq_Fct.R")
(varreport <- variantReport(allvar, vcf))[1:4,]
write.table(varreport, "results/varreport.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 29: Rvarseq.Rnw:690-692
###################################################
coding <- predictCoding(vcf, txdb, seqSource=fa)
coding[1:3,c(12,16:17)]


###################################################
### code chunk number 30: Rvarseq.Rnw:697-700
###################################################
source("Rvarseq_Fct.R")
(codereport <- codingReport(coding, txdb))[1:3,]
write.table(codereport, "results/codereport.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 31: Rvarseq.Rnw:712-715
###################################################
fullreport <- cbind(varreport, codereport[rownames(varreport),-1])
write.table(fullreport, "results/fullreport.xls", row.names=FALSE, quote=FALSE, sep="\t", na="")
fullreport[c(1,18),]


###################################################
### code chunk number 32: Rvarseq.Rnw:728-736
###################################################
library(VariantTools)
vr <- as(vcf, "VRanges")
varid <- paste(as.character(seqnames(vr)), ":", start(vr), "_", ref(vr), "/", alt(vr), sep="")
vrdf <- data.frame(row.names=varid, as.data.frame(vr))
vrdf <- vrdf[,c("totalDepth", "refDepth", "altDepth", "n.read.pos", "QUAL", "mean.quality")]
fullreport <- cbind(VARID=fullreport[,1], vrdf[rownames(fullreport),], fullreport[,-1])
fullreport[c(1,18),c(1:8,14)]
write.table(fullreport, "results/fullreport.xls", row.names=FALSE, quote=FALSE, sep="\t", na="")


###################################################
### code chunk number 33: Rvarseq.Rnw:761-770 (eval = FALSE)
###################################################
## library(SRAdb)
## startIGV("lm")
## sock <- IGVsocket()
## session <- IGVsession(files=c("results/SRR064154.fastq.bam", 
##                              "results/varianttools.vcf.bgz"), 
##                              sessionFile="session.xml", 
##                              genome="A. thaliana (TAIR10)")
## IGVload(sock, session)
## IGVgoto(sock, 'Chr5:6455')


###################################################
### code chunk number 34: Rvarseq.Rnw:783-790
###################################################
library(ggbio); library(GenomicAlignments)
ga <- readGAlignmentsFromBam(path(bfl[[1]]), use.names=TRUE, param=ScanBamParam(which=GRanges("Chr5", IRanges(4000, 8000))))
p1 <- autoplot(ga, geom = "rect")
p2 <- autoplot(ga, geom = "line", stat = "coverage")
p3 <- autoplot(vcf[seqnames(vcf)=="Chr5"], type = "fixed") + xlim(4000, 8000) + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank())
p4 <- autoplot(txdb, which=GRanges("Chr5", IRanges(4000, 8000)), names.expr = "gene_id")
tracks(Reads=p1, Coverage=p2, Variant=p3, Transcripts=p4, heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")


###################################################
### code chunk number 35: Exercise 2 (eval = FALSE)
###################################################
## bfl <- BamFileList(paste("./results/", as.character(targets[,1]), ".bam", sep=""), index=character())
## bfl
## ## Call variants, annotate them and store results for each sample in its own list component
## vrl <- sapply(basename(names(bfl)), function(x) NULL, simplify=FALSE)
## for(i in seq(along=bfl)) {
##     var <- callVariants(bfl[[i]], tally.param)
##     var <- var[totalDepth(var) == altDepth(var) & totalDepth(var)>=5 & values(var)$n.read.pos >= 5] 
##     sampleNames(var) <- "bwa"
##     vcf <- asVCF(var); rd <- rowData(vcf)
##     rownames(vcf) <- paste(as.character(seqnames(rd)), ":", start(rd), "_", values(rd)$REF, "/", values(rd)$ALT, sep="")
##     allvar <- locateVariants(rowData(vcf), txdb, AllVariants())
##     source("Rvarseq_Fct.R")
##     varreport <- variantReport(allvar, vcf)
##     writeVcf(vcf, "./results/variants.vcf")
##     vcf <- readVcf(file="results/variants.vcf", genome="ATH1")
##     coding <- predictCoding(vcf, txdb, seqSource=fa)
##     codereport <- codingReport(coding, txdb)
##     vrl[[i]] <- cbind(varreport, codereport[rownames(varreport),-1])[,-10]
##     colnames(vrl[[i]]) <- paste(colnames(vrl[[i]]), gsub("\\..*", "", names(vrl)[i]), sep="_")
##     print(i)
##     }
## ## Assemble results in a single matrix
## rowlabels <- unique(unlist(sapply(vrl, rownames, simplify=FALSE)))
## collabels <- unlist(sapply(vrl, colnames, simplify=FALSE))
## df <- as.data.frame(matrix(data=NA, nrow=length(rowlabels), ncol=length(collabels), dimnames=list(rowlabels, collabels)))
## for(i in seq(along=vrl)) {
##     df[rownames(vrl[[i]]), colnames(vrl[[i]])] <- vrl[[i]]
##     }
## df[1:20, 1:20]
## table(rowSums(!is.na(df[, grep("VARID", colnames(df))]))) # Variant occurrence in samples
## write.table(df, "results/fullreport4.xls", col.names=NA, quote=FALSE, sep="\t", na="")


###################################################
### code chunk number 36: Rvarseq.Rnw:845-846
###################################################
sessionInfo()


