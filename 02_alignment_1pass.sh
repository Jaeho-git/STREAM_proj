#!/bin/bash
#Alignment
Index_DIR=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/star_index
DATA=/home/jaeho/STREAM/RNA/00_Data/00_raw.data
Tools=/home/jaeho/Tools/STAR-2.7.10a/bin/Linux_x86_64

#source ./00_source
mkdir /home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass1

for i in $(ls $DATA | awk '$0 ~ ".fastq.gz"' | sed 's/_[12].fastq.gz//')
do
	mkdir /home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass1/$i
	${Tools}/STAR --genomeDir $Index_DIR \
--readFilesIn $DATA/${i}_1.fastq.gz $DATA/${i}_2.fastq.gz \
--runThreadN 12 \
--readFilesCommand zcat \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMtype None \
--outSAMmode None \
--outFileNamePrefix /home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass1/$i/$i. > /home/jaeho/STREAM/RNA/00_Data/01_processed.data/log/log_02_Alignment_pass1_$i
done

