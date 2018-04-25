#!/bin/bash

PREFIX="1yr"

source init.sh

for ID in "$@"
do
echo "<$ID>"
bamfile=bam/${PREFIX}/${ID}.sorted.bam
if [ -f "$bamfile" ] && [ -s "$bamfile" ]; then
    echo "exist non-empty $bamfile, skip"
    continue
fi
fastq-dump -O fastq ${ID}
./lib/STAR/source/STAR --genomeDir index_dir/index_22/  --runThreadN 8 --readFilesIn fastq/${ID}.fastq --outFileNamePrefix sam/${PREFIX}/${ID}
rm fastq/${ID}.fastq
mkdir bam/${PREFIX}
samtools view -S -b sam/${PREFIX}/${ID}Aligned.out.sam > bam/${PREFIX}/${ID}.raw.bam
samtools sort bam/${PREFIX}/${ID}.raw.bam -o bam/${PREFIX}/${ID}.sorted.bam
done

