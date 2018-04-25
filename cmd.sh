#!/bin/bash

# gen ref index
#nohup ./STAR/source/STAR  --runMode genomeGenerate --runThreadN 15 --genomeDir index_dir/ --genomeFastaFiles index_dir/GRCh38.p10.genome.fa > tmp.log &

nohup ./lib/STAR/source/STAR  --sjdbGTFfile index_dir/index_22/gencode.v28.annotation.gtf --runMode genomeGenerate --runThreadN 8 --genomeDir index_dir/index_22 --genomeFastaFiles index_dir/index_22/chr22.fa > nohup.log &

# download sra
#python lib/geofetch/geofetch/geofetch.py -i GSE81547 -m meta/

# get first 10 sample from first people
#find ncbi/public/sra/ -type f -name '*.sra' | sort | head -n 100000 | xargs -n 1 basename | cut -d '.' -f 1 > srr_list

# convert sra into fastq
#cat srr_list | xargs fastq-dump -O fastq

