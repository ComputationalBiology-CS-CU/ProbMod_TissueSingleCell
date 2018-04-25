#!/bin/bash
nohup sh -c "cat srr_list/1yr.txt | xargs sh run_align.sh" > nohup1.log &


