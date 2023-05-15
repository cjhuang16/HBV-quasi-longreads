#! /bin/sh
# filtering ccs reads carrying barcode sequence
# Usage: ./filter_out_barcode_carrying_reads.sc <SAMPLE>

#  May.2023
# Dr. Chih-Jen Huang
# Genomic Research Center, Academia Sinica, Taipei, Taiwan

SAMPLE="1810095"
# SAMPLE=$1

InD="./example"
OutD1="./CCS3/rmReadsWithBarcode/noBC/${SAMPLE}"
OutD2="./CCS3/rmReadsWithBarcode/innerBC_reads/${SAMPLE}"
Log="./CCS3/rmReadsWithBarcode/log"
barcode="./example/Barcode.seq.txt"


###################################
echo -e "                  \033[32mProcessing: $SAMPLE\033[0m             "
mkdir -p $OutD1 $OutD2 $Log

# extract reads name
zcat $InD/${SAMPLE}.ccs3.hifi_reads.fastq.gz | zgrep "^@m5" > ${Log}/${SAMPLE}.all.read.txt 2>>${Log}/${SAMPLE}.log.txt

# reads carrying barcode sequence
zcat $InD/${SAMPLE}.ccs3.hifi_reads.fastq.gz | zgrep -B1 -f ${barcode} | grep "^@m5" > ${Log}/${SAMPLE}.read.with.bc.txt 2>>${Log}/${SAMPLE}.log.txt

# reads without barcode
grep -v -f ${Log}/${SAMPLE}.read.with.bc.txt ${Log}/${SAMPLE}.all.read.txt > ${Log}/${SAMPLE}.noBC.txt 2>>${Log}/${SAMPLE}.log.txt

# output reads containing no barcode sequence
zcat $InD/${SAMPLE}.ccs3.hifi_reads.fastq.gz |zgrep -A3 -f ${Log}/${SAMPLE}.noBC.txt | sed '/^--$/d' |  gzip -c -  > $OutD1/${SAMPLE}.ccs3.hifi_reads.noBC.fastq.gz 2>>${Log}/${SAMPLE}.log.txt

# output reads with barcode sequence inside
zcat $InD/${SAMPLE}.ccs3.hifi_reads.fastq.gz |zgrep -A3 -f ${Log}/${SAMPLE}.read.with.bc.txt | sed '/^--$/d' |gzip -c -  > $OutD2/${SAMPLE}.ccs3.hifi_reads.innerBC.fastq.gz 2>>${Log}/${SAMPLE}.log.txt


