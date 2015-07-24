#!/bin/bash

module load bwa/0.7.10
module load gatk/3.3-0
module load abra/0.94
module load picard/1.121
module load samtools/1.1
module load varscan/2.3.7
module load sambamba
module load samblaster

patient='CCD1_53'
sample='tumor' # blood or tumor
read1='CCD1-53-Total-WES-1_S5_1.fastq.gz'
read2='CCD1-53-Total-WES-1_S5_2.fastq.gz'

## patient file directories
sample_dir="/data/s3/averafastq/patients/$patient"
alignment_dir="/data/storage/adam/patients/$patient/bam"
variant_dir="/data/storage/adam/patients/$patient/varscan_vcf"

mkdir -p /data/storage/adam/patients/$patient
mkdir -p /data/storage/adam/patients/$patient/varscan_vcf

## reference files
bwa_index="/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
bed="/data/database/GATK/beds/truseq_exome_targeted_regions.hg19.bed"
hg19_reference="/data/database/GATK/hg19/ucsc.hg19.fasta"

## abra realigner
#java -Xmx16G -jar `which abra.jar` \
#--in $bam \
#--out $out_dir/${name}\_realigned.bam \
#--ref $ref \
#--targets $bed \
#--threads 8 \
#--working abra_temp_dir > abra.log 2>&1

## picard mark duplicates
#java -Xmx16g -Djava.io.tmpdir=/tmp -jar `which MarkDuplicates.jar` \
#I=$alignment_dir/$patient\_$sample\_sorted.bam \
#O=$alignment_dir/$patient\_$sample\_WES_dedup.bam \
#M=$alignment_dir/$patient\_$sample\_WES_metrics.txt \
#ASSUME_SORTED=true \
#VALIDATION_STRINGENCY=LENIENT \

## samtools mpileup
#samtools mpileup -f $ref $alignment_dir/$patient\_$sample\_WES_dedup.bam > $alignment_dir/$patient\_$sample\_WES_dedup.pileup

## varscan2 call variants
VarScan.v2.3.7.jar pileup2snp $alignment_dir/$patient\_$sample\_WES_dedup.pileup --p-value 0.05 > $alignment_dir/$patient\_$sample\_WES_varscan.vcf

exit 0
