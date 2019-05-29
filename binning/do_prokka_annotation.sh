#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate bins with prokka
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/bins_of_interest/bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
cd /scratch/sogin/tmp.$JOB_ID/arc/
for i in *fa; do 
	prokka  --cpus 24 --kingdom Archaea --outdir "$i"_out $i; 
done
rm *fa
#
cd /scratch/sogin/tmp.$JOB_ID/bac/
for i in *fa; do 
	prokka --cpus 24 --outdir prokka_"$i"_result $i; 
done
rm *fa
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/bins_of_interest/prokka/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date