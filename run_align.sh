#!/bin/bash

PREFIX="1yr"

source init.sh
#./lib/STAR/source/STAR --genomeDir index_dir/index_dir/  --runThreadN 16 --readFilesIn fastq/${ID}.fastq --outFileNamePrefix sam/${ID}

for ID in "$@"
do
echo "<$ID>"
bamfile=bam/${PREFIX}/${ID}.sorted.bam
if [ -f "$bamfile" ] && [ -s "$bamfile" ]; then
    echo "exist non-empty $bamfile, skip"
    continue
fi
fastq-dump -O fastq ${ID}
#./lib/STAR/source/STAR --genomeDir index_dir/index_dir/  --runThreadN 16 --readFilesIn fastq/${ID}.fastq --outFileNamePrefix sam/${ID}
mkdir sam/${PREFIX}
./lib/STAR/source/STAR --genomeDir index_dir/index_22/  --runThreadN 8 --readFilesIn fastq/${ID}.fastq --outFileNamePrefix sam/${PREFIX}/${ID}
rm fastq/${ID}.fastq
mkdir bam/${PREFIX}
samtools view -S -b sam/${PREFIX}/${ID}Aligned.out.sam > bam/${PREFIX}/${ID}.raw.bam
samtools sort bam/${PREFIX}/${ID}.raw.bam -o bam/${PREFIX}/${ID}.sorted.bam
#samtools mpileup -uf index_dir/GRCh38.p10.genome.fa ${ID}.sorted.bam | bcftools call -mv -Ob -o ${ID}.bcf
#bcftools view ${ID}.bcf > ${ID}.vcf
done

IDs=""
IDfastq=""
for ID in "$@"
do
 IDfastq="$IDfastq fastq/$ID.fastq"    
 IDs="$IDs bam/$ID.sorted.bam" #//$ should be removed which is prefixed before var. Blank space before and after equal sign should be removed to run this code.   
done
echo $IDfastq
echo $IDs

#merge aligned BAM
#samtools merge bam/21yr.bam ${IDs}

#./lib/STAR/source/STAR --genomeDir index_dir/index_dir/  --runThreadN 16 --readFilesIn $IDfastq --outFileNamePrefix sam/21yr_
#samtools mpileup -uf index_dir/index_dir/GRCh38.p10.genome.fa ${IDs} | bcftools call -mv -Ob -o bcf/${ID}.bcf
#bcftools view bcf/${ID}.bcf > vcf/${ID}.vcf

