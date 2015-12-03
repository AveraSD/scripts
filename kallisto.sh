#!/bin/bash

index=/data/database/kallisto/gencode.v19.pc_transcripts_ercc92.idx
out=/data/storage/GeparQuinto/kallisto
cpu=$(grep -c "processor" /proc/cpuinfo)

module load kallisto

for i in $(find /data/storage/GeparQuinto/cleanfastq/ -type f -maxdepth 1 | sed 's/_..fastq.gz//' | sort | uniq)
do
  name=basename $i
  mkdir -p $out/$name
  kallisto quant -i $index -o $out/$name ${i}_1.fastq.gz ${i}_2.fastq.gz -b 100 -t $cpu 
done
