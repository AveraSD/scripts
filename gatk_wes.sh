#!/bin/bash


module load bwa/0.7.10
module load picard/1.121
module load samtools/1.1
module load gatk/3.3-0
module load R/3.1.1
module load sambamba
module load samblaster

patient='CCD1_53'
sample='tumor' # blood or tumor

mkdir -p /data/storage/adam/patients/$patient
mkdir -p /data/storage/adam/patients/$patient/bam
mkdir -p /data/storage/adam/patients/$patient/vcf
mkdir -p /data/storage/adam/patients/$patient/recal

## patient file directories
sample_dir="/data/s3/averafastq/patients/$patient"
alignment_dir="/data/storage/adam/patients/$patient/bam"
variant_dir="/data/storage/adam/patients/$patient/vcf"
recal_dir="/data/storage/adam/patients/$patient/recal"

## reference files
bwa_index="/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
RGR="@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:$patient\_$sample"
bed="/data/database/GATK/beds/truseq_exome_targeted_regions.hg19.bed"
hg19_reference="/data/database/GATK/hg19/ucsc.hg19.fasta"
dbsnp="/data/database/GATK/hg19/dbsnp_137.hg19.vcf"
mills="/data/database/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf"
g1000="/data/database/GATK/hg19/1000G_phase1.indels.hg19.vcf"
hapmap="/data/database/GATK/hg19/hapmap_3.3.hg19.vcf"
omni="/data/database/GATK/hg19/1000G_omni2.5.hg19.vcf"
KG="/data/database/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.vcf"

## create read group
java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which AddOrReplaceReadGroups.jar` \
I=$alignment_dir/$patient\_$sample\.sam \
O=$alignment_dir/$patient\_$sample\_sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=5000000 \
ID=1 \
LB=SRR \
PL=ILLUMINA \
PU=1234 \
SM=$sample

## picard mark duplicates
java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which MarkDuplicates.jar` \
I=$alignment_dir/$patient\_$sample\_sorted.bam \
O=$alignment_dir/$patient\_$sample\_dedup.bam \
M=$alignment_dir/$patient\_$sample\_metrics.txt \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT \

## index the bam file
java -Xmx16g -jar `which BuildBamIndex.jar` \
I=$alignment_dir/$patient\_$sample\_dedup.bam

### GATK
java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which GenomeAnalysisTK.jar` \
-T RealignerTargetCreator -nt 8 \
-R $hg19_reference \
-I $alignment_dir/$patient\_$sample\_dedup.bam \
-o $alignment_dir/$patient\_$sample\.intervals \
-L $bed \
-ip 100 \
-known $mills \
-known $g1000

#java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which GenomeAnalysisTK.jar` \
-T IndelRealigner \
-R $hg19_reference \
-I $alignment_dir/$patient\_$sample\_dedup.bam \
-targetIntervals $alignment_dir/$patient\_$sample\.intervals \
-o $alignment_dir/$patient\_$sample\_realigned.bam \
-known $mills \
-known $g1000

rm $alignment_dir/$patient\_$sample\_dedup.bam
rm $alignment_dir/$patient\_$sample\_dedup.bam.bai

#java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which GenomeAnalysisTK.jar` \
-T BaseRecalibrator -nct 8 \
-R $hg19_reference \
-I $alignment_dir/$patient\_$sample\_realigned.bam  \
-o $alignment_dir/$patient\_$sample\_recal_data.table \
-L $bed \
-ip 100 \
-knownSites $dbsnp \
-knownSites $mills \
-knownSites $g1000

java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which GenomeAnalysisTK.jar` \
-T PrintReads -nct 8 \
-R $hg19_reference \
-I $alignment_dir/$patient\_$sample\_realigned.bam \
-BQSR $alignment_dir/$patient\_$sample\_recal_data.table \
-o $alignment_dir/$patient\_$sample\_recal.bam

rm $alignment_dir/$patient\_$sample\.intervals
rm $alignment_dir/$patient\_$sample\_realigned.bam
rm $alignment_dir/$patient\_$sample\_realigned.bai
rm $alignment_dir/$patient\_$sample\_recal_data.table


## picard mark duplicates
java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which MarkDuplicates.jar` \
I=$alignment_dir/$patient\_$sample\_recal.bam \
O=$alignment_dir/dedup_$patient\_$sample\.bam \
M=$alignment_dir/$patient\_$sample\_metrics.txt \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT \

## index the bam file
java -Xmx16g -jar `which BuildBamIndex.jar` \
I=$alignment_dir/dedup_$patient\_$sample\.bam

## call variants
java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
-T HaplotypeCaller \
-nct 8 \
-R $hg19_reference \
--dbsnp $dbsnp \
-I $alignment_dir/dedup_$patient\_$sample\.bam \
-L $bed \
-ip 100 \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
-o $variant_dir/$patient\_$sample\.vcf

## recalibrate variant quality scores
java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
-T VariantRecalibrator \
-R $hg19_reference \
-input $variant_dir/$patient\_$sample\.vcf \
-nt 8 \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $KG \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile $recal_dir/$patient\_$sample\_recalibrate_SNP.recal \
-tranchesFile $recal_dir/$patient\_$sample\_recalibrate_SNP.tranches \
-rscriptFile $recal_dir/$patient\_$sample\_recalibrate_SNP_plots.R

java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
-T ApplyRecalibration \
-R $hg19_reference \
-input $variant_dir/$patient\_$sample\.vcf \
-tranchesFile $recal_dir/$patient\_$sample\_recalibrate_SNP.tranches \
-recalFile $recal_dir/$patient\_$sample\_recalibrate_SNP.recal \
-o $recal_dir/$patient\_$sample\_recalibrate_SNP.vcf \
--ts_filter_level 99.5 \
-mode SNP

#java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
#-T VariantRecalibrator \
#-R $hg19_reference \
#-input $recal_dir/$patient\_$sample\_recalibrate_SNP.vcf \
#-resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
#-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
#--maxGaussians 4 \
#--minNumBadVariants 5000 \
#-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
#-mode INDEL \
#-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
#-recalFile $recal_dir/$patient\_$sample\_recalibrate_INDEL.recal \
#-tranchesFile $recal_dir/$patient\_$sample\_recalibrate_INDEL.tranches \
#-rscriptFile $recal_dir/$patient\_$sample\_recalibrate_INDEL_plots.R
#
#java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
#-T ApplyRecalibration \
#-R $hg19_reference \
#-input $recal_dir/$patient\_$sample\_recalibrate_SNP.vcf \
#-tranchesFile $recal_dir/$patient\_$sample\_recalibrate_INDEL.tranches \
#-recalFile $recal_dir/$patient\_$sample\_recalibrate_INDEL.recal \
#-o $recal_dir/$patient\_$sample\_final.vcf \
#--ts_filter_level 99.0 \
#-mode INDEL

exit 0


