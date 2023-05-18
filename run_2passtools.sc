#!/bin/bash
## merge all bam files from the 1st run minimap2 alignment
## Usage: run_2passtools.sc <SAMPLE>
## Prerequisite software:
## 2passtools: https://github.com/bartongroup/2passtools
## samtools: http://www.htslib.org/
#
## Dr. Chih-Jen Huang
## Genomic Research Center, Academia Sinica, Taipei, Taiwan
## May.2023
#

ref="./example/HBV_isolate_SH1212-C5_complete_genome.fasta"
OutD="./2passtools"
InD="./minimap2_mapping/SH1212-C5_ref_splice-hq" 

mkdir -p $OutD


echo -e "\033[32m       >>> merging bam files\033[0m             "
# create a file list for samtools merge
ls ${InD}/*/*_minimap2_filter_final.clip10.sorted.bam |sed ':a;N;$!ba;s/\n/ /g' > ${OutD}/bam_list.txt

# merge bam file with samtools
# prepare a reading group file as rg.txt
samtools merge ${OutD}/merged.bam `cat ${OutD}/bam_list.txt`

samtools sort ${OutD}/merged.bam > ${OutD}/merged.sorted.bam
samtools index ${OutD}/merged.sorted.bam

echo -e "\033[32m       >>> ...finished\n           ...edit the header of the bam file further if you need or samtools will use the header from the first bam file\033[0m             "

# run 2passtools of 1st pass merged bam 
echo -e "\033[32m       >>> running 2passtools\033[0m             "

2passtools score -j 4  -d 10 -m 'GTAG|GCAG' -w 128 -k 6 -lt 0.1 -ht 0.6 --stranded -s 12345 -f $ref -o ${OutD}/merged_bam.2passtools.score.bed ${OutD}/merged.sorted.bam

# filter GTAG|GCAG plus strand with read counts > 100
awk '$4~/GTAG/ && $6~/+/ && $5>100' ${OutD}/merged_bam.2passtools.score.bed | cut -f-6 > ${OutD}/merged_bam.2passtools.score.plus.count100.bed 

awk '$4~/GCAG/ && $6~/+/ && $5>100' ${OutD}/merged_bam.2passtools.score.bed | cut -f-6 >> ${OutD}/merged_bam.2passtools.score.plus.count100.bed

# sort
sort -k2,2n -k3,3n ${OutD}/merged_bam.2passtools.score.plus.count100.bed > temp.txt && mv temp.txt ${OutD}/merged_bam.2passtools.score.plus.count100.bed

echo -e "\033[32m       >>> ...finished\n           ...output file: ${OutD}/merged_bam.2passtools.score.plus.count100.bed\033[0m             "

