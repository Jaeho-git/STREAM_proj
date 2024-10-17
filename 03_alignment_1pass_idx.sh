#!/bin/bash
#Indexing after 1pass alignment
DATA=/home/jaeho/STREAM/RNA/00_Data/00_raw.data
FASTA=/home/jaeho/Data/Reference/UCSC_resource/hg38.fa
Tools=/home/jaeho/Tools/STAR-2.7.10a/bin/Linux_x86_64

#source ./00_source
mkdir /home/jaeho/STREAM/RNA/00_Data/01_processed.data/star_aligned_pass1_index

Index_DIR=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/star_aligned_pass1_index
INPUT=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass1

for i in $(ls $DATA | awk '$0 ~ ".fastq.gz"' | sed 's/_[12].fastq.gz//' | sort -u)
do
    mkdir $Index_DIR/$i
	OUTPUT=$Index_DIR/$i
	${Tools}/STAR --runMode genomeGenerate --genomeDir $OUTPUT \
--genomeFastaFiles $FASTA \
--sjdbOverhang 100 \
--runThreadN 12 \
--readFilesCommand zcat \
--sjdbFileChrStartEnd $INPUT/$i/$i.SJ.out.tab
done
