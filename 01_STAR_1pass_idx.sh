#!/bin/bash
#source ./00_source
mkdir /home/jaeho/STREAM/RNA/00_Data/01_processed.data/log
mkdir /home/jaeho/STREAM/RNA/00_Data/01_processed.data/star_index
# STAR genome index : 1pass
Index_DIR=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/star_index
Ref_genome=/home/jaeho/Data/Reference/UCSC_resource/hg38.fa
Annotation=/home/jaeho/Data/Reference/Annotation/gencode.v46.annotation.gtf
Tools=/home/jaeho/Tools/STAR-2.7.10a/bin/Linux_x86_64

${Tools}/STAR --runMode genomeGenerate --genomeDir $Index_DIR --genomeFastaFiles $Ref_genome --sjdbOverhang 100 --sjdbGTFfile $Annotation --runThreadN 12 > /home/jaeho/STREAM/RNA/00_Data/01_processed.data/log/log_01_idx