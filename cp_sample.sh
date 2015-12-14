#!/bin/bash

sample=$1
path=/data/basespace/Projects/$sample/Samples
dest=/data/s3/averafastq/patients/$sample/

ss=()
counter=0
for i in $(find $path -maxdepth 1 -mindepth 1 -type d)
do
        ss[$counter]=$(basename $i | sed 's/_Rep.//gI')
        counter=$[$counter +1]
done

for i in $(echo "${ss[@]}" | tr '[:lower:]' '[:upper:]' | tr ' ' '\n' | sort -u | tr '\n' ' ')
do
        if [ $(echo $i | grep -cim1 DNA) -ge 1 ]; then
                if [ ! -f "${dest}/dna/${NAME1}.fastq.gz" ]; then
                        mkdir -p ${dest}/dna
                        cat $path/${i}_*/Files/*R1* > $dest/dna/${i}_1.fastq.gz &
                        cat $path/${i}_*/Files/*R2* > $dest/dna/${i}_2.fastq.gz
                else
                        echo "File exists"
                fi
        elif [ $(echo $i | grep -cim1 RNA) -ge 1 ]; then
                if [ ! -f "${dest}/rna/${NAME1}.fastq.gz" ]; then
                        mkdir -p ${dest}/rna
                        cat $path/${i}_*/Files/*R1* > $dest/rna/${i}_1.fastq.gz &
                        cat $path/${i}_*/Files/*R2* > $dest/rna/${i}_2.fastq.gz
                else
                        echo "File exists"
                fi
        else
                echo "Neither DNA or RNA in file name, don't know where to copy."
        fi
done

exit 0
