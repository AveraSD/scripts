#!/bin/bash

module load bwa/0.7.10
module load samtools/1.1
module load sambamba
module load samblaster
cpu=$(grep -c "processor" /proc/cpuinfo)

patient=$patient
sample=$sample # blood or tumor
read1=$read1
read2=$read2

## patient file directories
sample_dir="/data/s3/averafastq/patients/$patient"
alignment_dir="/data/storage/adam/patients/$patient/bam"

mkdir -p /data/storage/adam/patients/$patient
mkdir -p /data/storage/adam/patients/$patient/bam

## reference files
bwa_index="/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
RGR="@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:$patient\_$sample"
hg19_reference="/data/database/GATK/hg19/ucsc.hg19.fasta"

## align to hg19
bwa mem -M -t $cpu -R $RGR $bwa_index $sample_dir/$read1 $sample_dir/$read2 \
| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 25G --tmpdir /tmp -o $alignment_dir/$patient\_$sample\_WES_sorted.bam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 25G --tmpdir /tmp -o $alignment_dir/$patient\_$sample\_WES_sorted.bam /dev/stdin) \
| samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 25G --tmpdir /tmp -o $alignment_dir/$patient\_$sample\_WES_sorted.bam /dev/stdin

## mark duplicates
sambamba markdup -t $cpu $alignment_dir/$patient\_$sample\_WES_sorted.bam $alignment_dir/$patient\_$sample\_WES_dedup.bam

## index the bam file
sambama index -t $cpu $alignment_dir/$patient\_$sample\_WES_sorted.bam
