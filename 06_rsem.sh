#!/bin/bash

DATA=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/aligned_pass2
RSEM_dir=/home/jaeho/Tools/RSEM-1.3.3
RSEM_Anno=/home/jaeho/STREAM/RNA/00_Data/01_processed.data/rsem_hg38

#source ./00_source

mkdir /home/jaeho/STREAM/RNA/01_rsem_output

for i in $(ls $DATA | awk '$0 ~ ".Aligned.toTranscriptome.out.bam"' | sed s/.Aligned.toTranscriptome.out.bam// | sort -u)
do
	$RSEM_dir/rsem-calculate-expression --bam -p 60 --paired-end $DATA/${i}.Aligned.toTranscriptome.out.bam $RSEM_Anno/hg38 /home/jaeho/STREAM/RNA/01_rsem_output/${i}
done
