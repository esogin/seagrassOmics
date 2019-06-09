#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate bins run dbcan 
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 

rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/target_bins/out_bins/run_dbcan/ /scratch/sogin/tmp.$JOB_ID/
bins=$(echo *.faa)
for i in $bins; do 
python run_dbcan.py $i protein --out_dir ${i%%.faa}_output --dia_cpu 24 --hmm_cpu 24 --hot_cpu 24;
done
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/target_bins/out_bins/run_dbcan/ 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
