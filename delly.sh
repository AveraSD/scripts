#!/bin/bash

module load delly

#exclusions=human.hg19.excl.tsv
outdir=/data/storage/14MS10038/delly
reference=/data/database/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
tumor=/data/s3/averapatients/14MS10038/wgs-tumor/pat001t-ready.bam
germline=/data/s3/averapatients/14MS10038/wgs-blood/pat001n-ready.bam

mkdir -p  $outdir

delly -t DEL -o $out_dir/del.vcf -g $reference $tumor $germline

exit 0
