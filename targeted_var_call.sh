#!/usr/bin/bash

snpeff=/home/tmeissner/snpEff/snpEff.jar
snpsift=/home/tmeissner/snpEff/SnpSift.jar

# use test sample; TL-24-JZYHE7YQ

sampleID="$1"
#sampleID=TL-24-JZYHE7YQ   #TL-24-UIUW5KJJ
dir=/mnt/ResearchData/MEM/Tempus/$sampleID
out_dir=/mnt/Research/TempusRNA
tumorVCFfb=$dir/DNA/$sampleID.soma.freebayes.vcf
tumorVCFpin=$dir/DNA/$sampleID.soma.pindel.vcf
bamRNA=$dir/RNA/${sampleID}_T_sorted.bam
refGenome=/mnt/Research/CCB/database/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_nochr.fa

# creat targets region from vcf
cat $tumorVCFfb | tail -n +174 | awk '{FS="\t";OFS="\t";print $1,$2-1,$2,$3, etc}' > /tmp/$sampleID.fb.bed
cat $tumorVCFpin | tail -n +57 | awk '{FS="\t";OFS="\t";print $1,$2-1,$2,$3, etc}' > /tmp/$sampleID.pin.bed
cat /tmp/$sampleID.fb.bed /tmp/$sampleID.pin.bed > /tmp/$sampleID.bed

#create targets file
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $tumorVCFfb | bgzip -c > /tmp/$sampleID.als.fb.tsv.gz && tabix -s1 -b2 -e2 /tmp/$sampleID.als.fb.tsv.gz
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $tumorVCFpin | bgzip -c > /tmp/$sampleID.als.pin.tsv.gz && tabix -s1 -b2 -e2 /tmp/$sampleID.als.pin.tsv.gz
zcat /tmp/$sampleID.als.fb.tsv.gz /tmp/$sampleID.als.pin.tsv.gz | sort -k1,1 -k2,2n | bgzip -c > /tmp/$sampleID.als.tsv.gz && tabix -s1 -b2 -e2 /tmp/$sampleID.als.tsv.gz

## save input annotation
java -jar $snpsift extractFields $tumorVCFfb "CHROM" "POS" "ANN[0].GENE" "ANN[0].HGVS_P" | tail -n +2 > /tmp/$sampleID.in1.anno.tsv
java -jar $snpsift extractFields $tumorVCFpin "CHROM" "POS" "ANN[0].GENE" "ANN[0].HGVS_P" | tail -n +2 > /tmp/$sampleID.in2.anno.tsv
cat /tmp/$sampleID.in1.anno.tsv /tmp/$sampleID.in2.anno.tsv | sort -k1,1n -k2,2n | awk '!seen[$0]++' | awk -v OFS='\t' '$3 && !$4{ $4="NA" }1' | column -t > /tmp/$sampleID.in.anno.tsv

#call variants using bcftools
bcftools mpileup \
    --redo-BAQ \
    --min-BQ 30 \
    --per-sample-mF \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --threads 8 \
    -d 1000 \
    -Ou \
    -R /tmp/$sampleID.bed \
    -f $refGenome $bamRNA | 
bcftools call \
    --threads 8 \
    -A \
    -C alleles \
    -T /tmp/$sampleID.als.tsv.gz \
    -m \
    -Ov - |
 bcftools +fill-tags -- -t AF,FORMAT/VAF > /tmp/$sampleID.RNA.vcf 

## annotate variants with snpeff
java -Xmx8G -jar $snpeff eff hg19 /tmp/$sampleID.RNA.vcf > /tmp/$sampleID.RNA.anno.vcf

#convert to table output
vcftools --vcf /tmp/$sampleID.RNA.anno.vcf --extract-FORMAT-info "VAF" --out /tmp/$sampleID.RNA.anno.vaf.txt 
java -jar $snpsift extractFields /tmp/$sampleID.RNA.anno.vcf "CHROM" "POS" "ANN[0].GENE" "ANN[0].HGVS_P" | awk -v OFS='\t' '$3 && !$4{ $4="NA" }1' | column -t  > /tmp/$sampleID.RNA.anno.sift.txt 
paste /tmp/$sampleID.RNA.anno.vaf.txt.VAF.FORMAT /tmp/$sampleID.RNA.anno.sift.txt | awk -v OFS='\t' '{ print $1, $2, $6, $7, $3}'  | tail -n +2 | sort -k1,1n -k2,2n | column -t > /tmp/$sampleID.RNA.tbl.tsv

# get raw depth data at variant locations; convert to 0-based bed style
zcat /tmp/$sampleID.als.tsv.gz | awk -v OFS='\t' '{ print $1, $2-1}' | sort -k1,1n -k2,2n  > /tmp/$sampleID.pos.txt
sambamba depth base \
    -q=30 \
    -F "mapping_quality > 0 and not (duplicate or failed_quality_control or unmapped or secondary_alignment)" \
    --min-coverage=0 \
    -L /tmp/$sampleID.pos.txt \
    $bamRNA | tail -n +2  | column -t > /tmp/$sampleID.RNA.bamstats.txt

cat /tmp/$sampleID.RNA.bamstats.txt | awk -v OFS='\t' '{ print $1, $2+1, $3, $4, $5, $6, $7, $8, $9, $10}' | column -t > /tmp/$sampleID.RNA.bamstats.coord.txt


# join everyting & copy final results
awk '{ key = $1 FS $2 };
       NR == FNR { t[key] = $0; next };
       key in t { print t[key]; next };
       1' /tmp/$sampleID.RNA.tbl.tsv /tmp/$sampleID.in.anno.tsv | awk -v OFS='\t' '$4 && !$5{ $5="NA" }1' | column -t > /tmp/$sampleID.merge1.tsv

awk '{ key = $1 FS $2 };
       NR == FNR { t[key] = $0; next };
       key in t { print t[key]; next };
       1' /tmp/$sampleID.RNA.bamstats.coord.txt /tmp/$sampleID.in.anno.tsv | 
        awk -v OFS='\t' '$2 && !$3{ $3="NA" }1' |
        awk -v OFS='\t' '$3 && !$4{ $4="NA" }1' |
        awk -v OFS='\t' '$4 && !$5{ $5="NA" }1' |
        awk -v OFS='\t' '$5 && !$6{ $6="NA" }1' |
        awk -v OFS='\t' '$6 && !$7{ $7="NA" }1' |
        awk -v OFS='\t' '$7 && !$8{ $8="NA" }1' |
        awk -v OFS='\t' '$8 && !$9{ $9="NA" }1' |
        awk -v OFS='\t' '$9 && !$10{ $10="NA" }1' | 
        column -t > /tmp/$sampleID.merge2.tsv

paste /tmp/$sampleID.merge1.tsv /tmp/$sampleID.merge2.tsv | column -s "$(printf '\t')" -t | awk -v OFS='\t' '{$6=$7=""; print $0}' | column -t > /tmp/$sampleID.RNA.final.tsv
cp /tmp/$sampleID.RNA.final.tsv $out_dir

#clean tmp
rm /tmp/$sampleID*

exit 0
