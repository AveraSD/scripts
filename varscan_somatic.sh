#!/bin/bash

module load samtools/1.1
module load varscan/2.3.7

patient=$patient
normal_bam=$normal_bam # blood or tumor
tumor_bam=$tumor_bam

mkdir -p /data/s3/averapatients/$patient/varscan_vcf

## patient file directories
variant_dir="/data/s3/averapatients/$patient/varscan_vcf"

## reference files
ref="/data/database/GATK/hg19/ucsc.hg19.fasta"

## samtools mpileup
samtools mpileup -f $ref $normal_bam > /tmp/${normal_bam##*/}.pileup
samtools mpileup -f $ref $tumor_bam > /tmp/${tumor_bam##*/}.pileup

## varscan2 call variants
VarScan.v2.3.7.jar somatic \
/tmp/${normal_bam##*/}.pileup \
/tmp/${tumor_bam##*/}.pileup \
--output-vcf 1 \
> $variant_dir/$patient\_varscan_somatic.vcf

#rm /tmp/*.pileup

exit 0
