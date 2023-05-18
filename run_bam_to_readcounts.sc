# run bam-readcount from sorted.bam
# Usage: run_bam_to_readcounts.sc <SAMPLE>
# Usage: run_2passtools.sc <SAMPLE>
## Prerequisite software:
# bam-readcount: https://github.com/genome/bam-readcount
#
## Dr. Chih-Jen Huang
## Genomic Research Center, Academia Sinica, Taipei, Taiwan
## May.2023

# exit when command failed
set -e

SAMPLE="1810095"
#SAMPLE=$1 # or use args

InD="./minimap2_mapping/SH1212-C5_ref_splice-hq-realign/${SAMPLE}"
OutD="./minimap2_mapping/SH1212-C5_ref_splice-hq-realign/${SAMPLE}/bam-readcount"
ref="./example/HBV_isolate_SH1212-C5_complete_genome.fasta"

echo -e "                  \033[32mProcessing: ${SAMPLE}\033[0m             "

mkdir -p $OutD

bamfile="_SH1212-C5_minimap2_filter_final.clip10.sorted.bam"

[ -f $InD/${SAMPLE}${bamfile} ] && bam-readcount -w1 -f $ref $InD/${SAMPLE}${bamfile} JX661488.1 > $OutD/${SAMPLE}.bam-readcount 2>$OutD/err.${SAMPLE}.bam-readcount || echo -e "\033[31mWarnning: bamfile ${SAMPLE}${bamfile} ...not found!\033[0m"
# JX661488.1: chromosome name of the reference genome
#
#
