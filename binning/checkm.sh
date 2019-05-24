#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
#
# run checkm
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/checkm/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
#checkm tree bins/ checkm_tree/ -t 24 -x .fa
#checkm tree_qa ./checkm_tree -o 2 -f tree_qa
#checkm lineage_set ./checkm_tree markers
#checkm analyze markers ./bins checkm_analyze -t 24 -x .fa
checkm qa markers ./checkm_analyze -o 2 -t 24 -f checkm_qa_results.txt
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/checkm/ 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
