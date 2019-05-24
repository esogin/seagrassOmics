#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 1
#$ -V
#$ -q main.q@@himem
#
# Assembly script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

./make_hit_count_tables.R

echo "job finished: "
date

