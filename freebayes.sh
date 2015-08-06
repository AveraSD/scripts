#!/bin/bash

module load samtools/1.1
module load vcflib
module load freebayes/0.9.21
cpu=$(grep -c "processor" /proc/cpuinfo)

patient='CCD1_53'
sample='tumor' # blood or tumor


## patient file directories
alignment_dir="/data/s3/averapatients/$patient/alignments"
variant_dir="/data/s3/averapatients/$patient/freebayes_vcf"

mkdir -p /data/s3/averapatients/$patient/freebayes_vcf

## reference files
hg19_reference="/data/database/GATK/hg19/ucsc.hg19.fasta"

## call variants with freebayes
freebayes \
    -f $hg19_reference \
    $alignment_dir/$patient\_$sample\_WES_realigned_sorted.bam \
     > $variant_dir/$patient\_$sample\_freebayes.vcf

exit 0
