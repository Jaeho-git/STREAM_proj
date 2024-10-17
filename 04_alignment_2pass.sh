#!/bin/bash
#Alignment

Index_DIR=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/star_aligned_pass1_index
DATA=/home/jaeho/STREAM/RNA/00_Data/00_raw.data
Annotation=/home/jaeho/Data/Reference/Annotation/gencode.v46.annotation.gtf
Tools=/home/jaeho/Tools/STAR-2.7.10a/bin/Linux_x86_64
#source ./00_source
mkdir /home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass2

for i in $(ls $DATA | awk '$0 ~ ".fastq.gz"' | sed 's/_[12].fastq.gz//' | sort -u)
do
	INPUT=$Index_DIR/$i
$Tools/STAR --genomeDir $INPUT \
--readFilesIn $DATA/${i}_1.fastq.gz $DATA/${i}_2.fastq.gz \
--runThreadN 10 \
--readFilesCommand zcat \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 10000000000 \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--sjdbGTFfile $Annotation \
--outSAMstrandField intronMotif \
--outSAMattributes NH HI NM MD AS XS \
--outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMheaderHD @HD VN:1.4 \
--outFileNamePrefix /home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass2/$i. > /home/jaeho/STREAM/RNA/00_Data/01_processed.data/log/log_04_Alignment_pass2_$i
done

