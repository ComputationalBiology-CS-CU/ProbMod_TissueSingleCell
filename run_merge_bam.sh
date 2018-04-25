#!/bin/bash

PREFIX="1yr"

source init.sh

IDs=""
for ID in "$@"
do
 IDs="$IDs bam/$PREFIX/$ID.sorted.bam"
done
echo $IDs

samtools merge -f bam/$PREFIX.merged.bam $IDs


