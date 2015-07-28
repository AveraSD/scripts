#!/bin/bash

read1=/data/s3/averafastq/everything_else/T47D-Dnase_S5_1.fastq.gz
read2=/data/s3/averafastq/everything_else/T47D-Dnase_S5_2.fastq.gz
genomeDir=/data/database/Homo_sapiens/UCSC/hg19/star_2.4.2a_genome
refgtf=/data/database/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
threads=32
rgid="1234"
rgsm="NA18238"
rgpl="ILLUMINA"
rglb="TrueSeq"
scratch=/tmp

STAR \
	--genomeDir $genomeDir \
	--readFilesIn $read1 $read2 \
	--readFilesCommand zcat \
	--twopassMode Basic \
	--quantMode TranscriptomeSAM GeneCounts \
	--sjdbGTFfile $refgtf \
	--alignIntronMax 200000 \
	--alignMatesGapMax 200000 \
	--outFilterMismatchNoverLmax 0.04 \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattrRGline ID:$rgid SM:$rgsm PL:$rgpl LB:$rglb \
	--outSAMstrandField intronMotif \
	--runThreadN $threads \
	--chimSegmentMin 25 \
	--chimJunctionOverhangMin 25 \
	--outTmpDir $scratch/$rgsm

rm -rf $scratch/$rgsm

exit 0
