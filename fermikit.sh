#!/bin/bash

module load fermikit/0.13

sample="14MS10038"
read1="/data/s3/averafastq/patients/14MS10038/tumor-r1.fq.gz"
read2="/data/s3/averafastq/patients/14MS10038/tumor-r2.fq.gz"
bwa_index="/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
out="/data/storage/fermikit"

cpu=$(grep -c "processor" /proc/cpuinfo)

cat $read1 $read2 > /tmp/$sample.fastq.gz

# assembly reads into unitigs (-s specifies the genome size and -l the read length)
fermi2.pl unitig -s3g -t$cpu -l75 -p $sample /tmp/$sample.fastq.gz > $out/$sample.mak
make -f $out/$sample.mak

# call small variants and structural variations
run-calling -t${cpu} $bwa_index $out/$sample.mag.gz | sh

# clean
rm /tmp/$sample.fastq.gz

exit 0
