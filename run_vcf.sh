#!/bin/bash

# input: prefix list
# e.g:
OUT_PREFIX="out"
mkdir $OUT_PREFIX
IDs=""
for ID in "$@"
do
 IDs="$IDs bam/$ID.merged.bam"
done
echo $IDs

samtools mpileup -uf ../index_dir/index_22/chr22.fa ${IDs} | bcftools call -mv -Ob -o out/${OUT_PREFIX}.bcf
bcftools view out/${OUT_PREFIX}.bcf > out/${OUT_PREFIX}.vcf

