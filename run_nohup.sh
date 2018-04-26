#!/bin/bash

array=( 1yr 5yr 6yr 21yr 22yr 38yr 44yr 54yr )
for i in "${array[@]}"
do
    echo $i
    nohup PREFIX=$i sh -c "cat srr_list/$i.txt | xargs sh run_align.sh" > nohup.log &
done



