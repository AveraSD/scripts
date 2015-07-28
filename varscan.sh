#!/bin/bash

module load samtools/1.1
module load varscan/2.3.7

patient='CCD1_53'
sample='tumor' # blood or tumor


## patient file directories
sample_dir="/data/s3/averafastq/patients/$patient"
alignment_dir="/data/s3/averapatients/$patient/bam"
variant_dir="/data/s3/averapatients/$patient/varscan_vcf"

mkdir -p /data/s3/averapatients/$patient/varscan_vcf

## reference files
ref="/data/database/GATK/hg19/ucsc.hg19.fasta"

## samtools mpileup
samtools mpileup -f $ref $alignment_dir/$patient\_$sample\_varscan_realigned.bam > $alignment_dir/$patient\_$sample\_varscan_realigned.pileup

## varscan2 call variants
VarScan.v2.3.7.jar pileup2snp $alignment_dir/$patient\_$sample\_varscan_realigned.pileup --p-value 0.05 > $variant_dir/$patient\_$sample\_varscan.vcf

exit 0
