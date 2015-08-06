#!/bin/bash

module load mutect/1.1.7

patient=$patient
normal_bam=$normal_bam
tumor_bam=$tumor_bam

mkdir -p /data/s3/averapatients/$patient
variant_dir=/data/s3/averapatients/$patient/gatk_vcf

hg19_reference="/data/database/GATK/hg19/ucsc.hg19.fasta"
dbsnp="/data/database/GATK/hg19/dbsnp_137.hg19.vcf"
cosimc="/data/database/cosmic/CosmicCodingMuts_v68_wchr_sort.vcf"
bed="/data/database/GATK/beds/truseq_exome_targeted_regions.hg19.bed"

java -Xmx2g -jar `which mutect-1.1.7.jar` \
--analysis_type MuTect \
--reference_sequence $hg19_reference \
--cosmic $cosmic \
--dbsnp $dbsnp \
--intervals $bed \
--input_file:normal $normal_bam \
--input_file:tumor $tumor_bam \
--vcf $variant_dir/$patient\_mutect.vcf \
-rf BadCigar \
--coverage_file $variant_dir/$patient\_coverage.wig.txt

exit 0
