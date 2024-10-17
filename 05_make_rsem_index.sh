#!/bin/bash

TOOLs=/Users/aimedbio-a10174/Tools
Ref_genome=/Users/aimedbio-a10174/Downloads/STREAM/Data/reference
Annotation=/Users/aimedbio-a10174/Downloads/STREAM/Data/annotation

mkdir rsem_hg38

$TOOLs/RSEM/rsem-prepare-reference -p 12 --star --star-path /usr/local/bin \
--gtf $Annotation/gencode.v46.annotation.gtf \
$Ref_genome/hg38.fa \
rsem_hg38/hg38
