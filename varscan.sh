#!/bin/bash

module load samtools/1.1
module load varscan/2.3.7

patient='CCD1_53'
sample='blood' # blood or tumor


## patient file directories
alignment_dir="/data/s3/averapatients/$patient/alignments"
variant_dir="/data/s3/averapatients/$patient/varscan_vcf"

mkdir -p /data/s3/averapatients/$patient/varscan_vcf

## reference files
ref="/data/database/GATK/hg19/ucsc.hg19.fasta"

## samtools mpileup
samtools mpileup -f $ref $alignment_dir/$patient\_$sample\_WES_realigned_sorted.bam > /tmp/$patient\_$sample\_WES_realigned_sorted.pileup

## varscan2 call variants
#VarScan.v2.3.7.jar pileup2snp /tmp/$patient\_$sample\_WES_realigned_sorted.pileup --p-value 0.05 > $variant_dir/$patient\_$sample\_varscan_snp.vcf
VarScan.v2.3.7.jar pileup2indel /tmp/$patient\_$sample\_WES_realigned_sorted.pileup --p-value 0.05 > $variant_dir/$patient\_$sample\_varscan_indel.vcf

rm /tmp/$patient\_$sample\_WES_realigned_sorted.pileup

exit 0
