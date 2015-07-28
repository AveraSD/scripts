#!/bin/bash

module load bwa/0.7.10
module load samtools/1.1
module load sambamba
module load samblaster
module load abra/0.94
cpu=$(grep -c "processor" /proc/cpuinfo)

patient=$patient
sample=$sample # blood or tumor
read1=$read1
read2=$read2
alignment_dir=$alignment_dir

mkdir -p /data/s3/averafastq/$patient
mkdir -p /data/s3/averapatients/$patient/bam
mkdir -p "$alignment_dir"

## reference files
bwa_index="/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
RGR="@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:$patient\_$sample"
hg19_reference="/data/database/GATK/hg19/ucsc.hg19.fasta"
bed="/data/database/GATK/beds/truseq_exome_targeted_regions.hg19.bed"

## align to hg19
bwa mem -M -t $cpu -R $RGR $bwa_index $read1 $read2 \
| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 25G --tmpdir /tmp -o $alignment_dir/$patient\_$sample\_WES_sorted.bam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 25G --tmpdir /tmp -o $alignment_dir/$patient\_$sample\_WES_sorted.bam /dev/stdin) \
| samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 25G --tmpdir /tmp -o $alignment_dir/$patient\_$sample\_WES_sorted.bam /dev/stdin

## mark duplicates
sambamba markdup -t $cpu $alignment_dir/$patient\_$sample\_WES_sorted.bam $alignment_dir/$patient\_$sample\_WES_dedup.bam

## index the bam file
sambamba index -t $cpu $alignment_dir/$patient\_$sample\_WES_dedup.bam

## realign with abra
java -Xmx4G -jar `which abra.jar` \
--in $alignment_dir/$patient\_$sample\_WES_dedup.bam \
--out $alignment_dir/$patient\_$sample\_WES_realigned.bam \
--ref $ref \
--targets $bed \
--threads 8 \
--working /tmp/abra_temp_dir > abra.log 2>&1
