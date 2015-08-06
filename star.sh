#!/bin/bash
module load star/2.4.2a

sample=$sample
read1=$read1
read2=$read2
alignment_dir=$alignment_dir
refgtf=$refgtf

mkdir -p "$alignment_dir"/$sample

cpu=$(grep -c "processor" /proc/cpuinfo)
genomeDir=/data/database/Homo_sapiens/UCSC/hg19/star_2.4.2a_genome

rgid="1234"
rgsm=$(echo $sample)
rgpl="ILLUMINA"
rglb="TrueSeq"
scratch=/tmp

# make sure temp is empty
if [ -d $scratch/$rgsm ]
  then
    rm -rf $scratch/$rgsm
fi

mkdir -p $scratch/$rgsm && cd $scratch/$rgsm

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
	--runThreadN $cpu \
	--chimSegmentMin 25 \
	--chimJunctionOverhangMin 25 \
	--outFileNamePrefix $sample.

cp $scratch/$rgsm/* $alignment_dir/$sample/

#clean up
rm -rf $scratch/$sample
rm -rf $scratch/$rgsm

exit 0
