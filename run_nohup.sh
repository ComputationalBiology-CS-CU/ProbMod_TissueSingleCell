#!/bin/bash
nohup sh -c "cat srr_list/22yr.txt | xargs sh run_align.sh" > nohup.log &


