#!/bin/bash

module load samtools/1.1
module load vcflib
module load freebayes/0.9.21
cpu=$(grep -c "processor" /proc/cpuinfo)

patient=$patient
normal_bam=$normal_bam # blood or tumor
tumor_bam=$tumor_bam

mkdir -p /data/s3/averapatients/$patient/freebayes_vcf

## patient file directories
variant_dir="/data/s3/averapatients/$patient/freebayes_vcf"

## reference files
ref="/data/database/GATK/hg19/ucsc.hg19.fasta"

## call variants with freebayes
freebayes \
-f $ref \
--pooled-continuous \
--pooled-discrete \
-F 0.03 \
-C 2 \
$tumor_bam \
$normal_bam \
>$variant_dir/$patient\_freebayes_somatic.vcf

exit 0
