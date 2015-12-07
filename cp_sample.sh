
#!/bin/bash

sample=$1

for i in $(find /data/basespace/Projects/$sample/Samples -maxdepth 1 -mindepth 1 -type d)
do
        NAME1=$(ls $i/Files | head -n 1 | sed 's/.fastq.gz//' | sed 's/_S.*/_1/')
        NAME2=$(ls $i/Files | head -n 1 | sed 's/.fastq.gz//' | sed 's/_S.*/_2/')

        if [ $(echo $NAME1 | grep -cim1 DNA) -ge 1 ]; then
                if [ ! -f "/data/s3/averafastq/patients/$sample/dna/${NAME1}.fastq.gz" ]; then
                        mkdir -p /data/s3/averafastq/patients/$sample/dna 
                        cat $i/Files/*R1* > /data/s3/averafastq/patients/$sample/dna/${NAME1}.fastq.gz &
                        cat $i/Files/*R2* > /data/s3/averafastq/patients/$sample/dna/${NAME2}.fastq.gz
                else
                        echo "File exists"
                fi
        elif [ $(echo $NAME1 | grep -cim1 RNA) -ge 1 ]; then
                if [ ! -f "/data/s3/averafastq/patients/$sample/rna/${NAME1}.fastq.gz" ]; then
                        mkdir -p /data/s3/averafastq/patients/$sample/rna 
                        cat $i/Files/*R1* > /data/s3/averafastq/patients/$sample/rna/${NAME1}.fastq.gz &
                        cat $i/Files/*R2* > /data/s3/averafastq/patients/$sample/rna/${NAME2}.fastq.gz
                else
                        echo "File exists"
                fi
        else
                echo "Neither DNA or RNA in file name, don't know where to copy."
        fi
done

exit 0
