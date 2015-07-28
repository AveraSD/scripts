#!/bin/bash

module load fermikit/0.13

sample="CCD1_53-4"
read1="/data/s3/averafastq/patients/CCD1_53/CCD1-53-Total-WES-4_S8_1.fastq.gz"
read2="/data/s3/averafastq/patients/CCD1_53/CCD1-53-Total-WES-4_S8_2.fastq.gz"
bwa_index="/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
out="/data/storage/fermi"

cpu=$(grep -c "processor" /proc/cpuinfo)

cat $read1 $read2 > /tmp/$sample.fastq.gz

# assembly reads into unitigs (-s specifies the genome size and -l the read length)
fermi2.pl unitig -s3g -t$cpu -l75 -p /tmp/$sample /tmp/$sample.fastq.gz > /tmp/$sample.mak
make -f /tmp/$sample.mak

# call small variants and structural variations
run-calling -t${cpu} $bwa_index /tmp/$sample.mag.gz | sh

# clean
rm /tmp/$sample.fastq.gz

# move to output
mv /tmp/$sample.* $out

exit 0
