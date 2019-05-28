#!/bin/bash
set -x
source $MODULESHOME/init/bash

module load bcftools

delly=/home/tobias/bin/delly_v0.8.1_linux_x86_64bit
hg19excl=/home/tobias/bin/human.hg19.excl.tsv
out=/data/storage/ikmb_gastric/delly
ref=/data/database/Homo_sapiens/b37decoy/hs37d5.fa

tbam=$1
gbam=$2
name=$(basename "$tbam" .bam)

tid=$(basename "$tbam" _realigned_sorted.bam)
gid=$(basename "$gbam" _realigned_sorted.bam)

awk -v tid="$tid" -v gid="$gid" 'BEGIN {print tid "\t" "tumor" "\n" gid "\t" "control"}' > $out/$tid.txt

$delly call \
	-x $hg19excl \
	-o $out/$name.bcf \
	-g $ref \
	-q 20 \
	-s 15 \
	$1 $2

$delly filter \
	-f somatic \
	-o $out/$name.somatic.bcf \
	-p \
	-s $out/$tid.txt \
	$out/$name.bcf


bcftools view $out/$name.somatic.bcf > $out/$name.somatic.vcf

exit 0
