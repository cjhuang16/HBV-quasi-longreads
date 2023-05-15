#!/bin/bash
#run minimap2 from file name list , when file name are too different
# ../../../sh/m1-minimap2_align_1st-round-SH1212-c5.sc <$1>

source /home/cjh/anaconda3/etc/profile.d/conda.sh
conda activate minimap2

ref="./example/HBV_isolate_SH1212-C5_complete_genome.fasta"
refNAME="SH1212-C5"
SAMPLE="1810095"
#SAMPLE=$1
OutD="./minimap2_mapping/${refNAME}_ref_splice-hq"
InD="./CCS3/rmReadsWithBarcode/noBC/${SAMPLE}"

RUNID="RUNID"

echo " ref= ${ref}"

mkdir -p $OutD

echo -e "\033[32mProcessing: $1\033[0m             "

mkdir -p $OutD/${SAMPLE}


# run command:
## minimap2 splice:hq mode
echo -e "\033[32m       >>> minimap2 mapping\033[0m             "
minimap2 -a --cs=long -x splice:hq -t 10 --secondary=no --end-bonus 5 -R "@RG\tID:${RUNID}\tPL:Pacbio_Sequel_ccs\tSM:${SAMPLE}\tRE:${refNAME}" ${ref} $InD/${SAMPLE}.ccs3.hifi_reads.noBC.fastq.gz -o $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.sam

# sam to bam
# -F0x900: exclude supplementary alignment (github.com/lh3/minimap2/blob/master/FAQ.md)
samtools view -F0x904 -Sb -h $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.bam
samtools view -f0x900 -Sb -h $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.f0x900.bam

#extract umapped reads
mkdir -p $OutD/${SAMPLE}/unmapped/
samtools view -f 4 -h $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.sam > $OutD/${SAMPLE}/unmapped/${SAMPLE}_${refNAME}_minimap2.unmapped.sam

# convert samfile
samtools view -h $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.bam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam

# sort bam files. The command is slightly different than the old samtools version
samtools sort $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.bam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sorted.bam

# make index of the sorted bam
samtools index $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sorted.bam

#########################
# Post alignment filter #
#########################

# Position and length
echo -e "\033[32m       >>> filtering: position and length\033[0m       "

#Chromosome name to filter for
CHROM="JX661488.1"

# keep mapped reads which start in between of $from and $to
from="50" # mapping results from 26  #2966 (type C control show non-specific amplification from here)
to="3200"  # mapping results end with 

[ -f $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam ] && echo -e "file: $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam  ...found" || (echo -e " \033[31m  file $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam not found. Please check\033[0m"; exit 1)

#filter by chromosome
cat $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam | awk -v CHROM=${CHROM} 'BEGIN{FS="\t";} {if($3==CHROM) print $0}' > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}.sam

#filter by given target (start and end)

cat $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}.sam | awk -v from=${from} -v to=${to} '{FS="\t"} $4<=from' from=$from to=$to > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam.left

cat $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}.sam | awk -v from=${from} -v to=${to} '{FS="\t"} $4>=to' from=$from to=$to > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam.right

cat $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam.left $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam.right > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam

#paste the original sam file header to the new sam file
head -n 100 $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam | grep  "^@" | cat - $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam > tmp.txt && mv tmp.txt $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam

#convert sam to bam
samtools view -bS $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.bam
# create sorted bam
samtools sort $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.bam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sorted.bam
#create bai
samtools index $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sorted.bam 

######################################################
# filter the tail through CIGAR tag
# calculate CIGAR length
samtools view $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sorted.bam | awk '{print $1","$6}' | sed 's/M/+/g' | sed 's/N/+/g' | sed 's/S/+/g'| sed 's/D/+/g'| sed 's/I/+/g'| sed 's/[0-9]*H/0+/g' | sed 's/+$//g' > temp.csv

for CIGARSUM in `awk -F, '{print $2}' temp.csv`; do echo $(($CIGARSUM)) >> temp2.csv; done

paste -d "," temp.csv temp2.csv | cut -d "," -f 1,3 > $OutD/${SAMPLE}/${SAMPLE}.cigar.lengthsum.csv

rm temp.csv temp2.csv
# CIGAR length less than 3000
awk -F, '{if ($2 <3000){print $1}}' $OutD/${SAMPLE}/${SAMPLE}.cigar.lengthsum.csv > $OutD/${SAMPLE}/${SAMPLE}.qname.cigarSum.less3000.txt

grep -v -f $OutD/${SAMPLE}/${SAMPLE}.qname.cigarSum.less3000.txt $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam

samtools view -bS $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.bam


# polyA tail
awk '{print $1"\t"$10}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -e 'TTCAAAAAAAAA' -e 'TGGAAAAAAAAA' -e 'AGCAAAAAAAAA' -e 'ATCAAAAAAAAA' -e 'TACAAAAAAAAA' | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.polyA.txt
# polyA longer than 12
awk '{print $1"\t"$10}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -e 'AAAAAAAAAAAA'  | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.polyA12.txt

awk '{if ($4 < 10){print $1}}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -v "^@" > $OutD/${SAMPLE}/${SAMPLE}.qname.mapStartBefore10.txt

awk '{print $1"\t"$10}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -e 'CCAGCACCATGCAACTTTTTCACCTCTGCCTAATCATCTC' | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.primerRFfusion.txt

# filter out reads with gap open within primer design region (alignment start from position 10: trimmed bacode and 6 primer
awk '{if ($4 == 10){print $1"\t"$6}}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -e $'\t'1M -e $'\t'2M -e $'\t'3M -e $'\t'4M -e $'\t'5M -e $'\t'6M -e $'\t'7M -e $'\t'8M -e $'\t'9M  -e $'\t'10M -e $'\t'11M -e $'\t'12M -e $'\t'13M -e $'\t'14M -e $'\t'15M -e $'\t'16M | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.FPrimerbreak.txt

# record read names with gap open at the very begining 
awk '{if ($4 == 10){print  $1"\t"$6}}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -e $'\t'17M -e $'\t'18M -e $'\t'19M -e $'\t'20M -e $'\t'21M -e $'\t'22M -e $'\t'23M -e $'\t'24M -e $'\t'25M -e $'\t'26M -e $'\t'27M -e $'\t'28M -e $'\t'29M -e $'\t'30M | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.start10.16-30M_not_removed.txt

awk '{print $1"\t"$6}' $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam | grep -e N1.M$ -e N2.M$ -e N1.M.S$ |cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.RPrimerbreak.txt

###########
cat $OutD/${SAMPLE}/${SAMPLE}.qname.polyA12.txt $OutD/${SAMPLE}/${SAMPLE}.qname.mapStartBefore10.txt $OutD/${SAMPLE}/${SAMPLE}.qname.FPrimerbreak.txt $OutD/${SAMPLE}/${SAMPLE}.qname.RPrimerbreak.txt $OutD/${SAMPLE}/${SAMPLE}.qname.primerRFfusion.txt | sort |uniq > $OutD/${SAMPLE}/${SAMPLE}.qname.to_be_remove.txt

grep -v -f $OutD/${SAMPLE}/${SAMPLE}.qname.to_be_remove.txt $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.sam

samtools view -bS $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.bam

samtools sort $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.bam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.sorted.bam

samtools index $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.sorted.bam

# clean intermediate file
#rm $samchromfile ${samchromtargetfile}.left ${samchromtargetfile}.right ${samchromtargetfile} ${bamchromtargetfile} ${sortedbamchromtargetfile} ${sortedbamchromtargetfile}.bai  ${samfilterlengthfile} ${SAMPLE}.cigar.lengthsum.csv

# filter out reads with clipping at head and tail

echo -e "\033[32m       >>> clipping filtering\033[0m             "

# find soft clip >= 10 at both end
samtools view $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.bam | awk '{print $1"\t"$6}' | grep -e $'\t'..S -e $'\t'...S -e $'\t'....S -e M....S$ -e M...S$ -e M..S$ | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.clip-s10.txt
# find hard clip >= 10 at both end
samtools view $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.bam | awk '{print $1"\t"$6}' | grep -e $'\t'..H -e $'\t'...H -e $'\t'....H -e M....H$ -e M...H$ -e M..H$ | cut -f1 > $OutD/${SAMPLE}/${SAMPLE}.qname.clip-h10.txt

cat $OutD/${SAMPLE}/${SAMPLE}.qname.clip-s10.txt $OutD/${SAMPLE}/${SAMPLE}.qname.clip-h10.txt > $OutD/${SAMPLE}/${SAMPLE}.qname.clip10.txt

grep -v -f $OutD/${SAMPLE}/${SAMPLE}.qname.clip10.txt $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.sam

samtools view -bS $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.sam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.bam

samtools sort $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.bam > $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.sorted.bam

samtools index $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.sorted.bam


# clean intermediate file
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.sam 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2.uniq.bam
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.sam
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}.sam 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam.left 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam.right 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sam 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.bam 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sorted.bam 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_${CHROM}_filter_${from}_to_${to}.sorted.bam.bai 
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.sam 
rm $OutD/${SAMPLE}/${SAMPLE}.cigar.lengthsum.csv
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_cigard.length_final.sam
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.bam
rm $OutD/${SAMPLE}/${SAMPLE}_${refNAME}_minimap2_filter_final.clip10.bam
echo -e "                  \033[32m          ... ${SAMPLE} run finished\033[0m    " 

