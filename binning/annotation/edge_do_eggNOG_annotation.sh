#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 12
#$ -V
#$ -q main.q
## Annotate bins with eggNOG
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/edge_bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/bins/
bins=$(echo *fa)

cd /scratch/sogin/tmp.$JOB_ID/

for i in $bins; 
	do
		cd /scratch/sogin/tmp.$JOB_ID/
		# run prodigial 
		prodigal -i bins/$i -o pr_run/genes/${i%%.fa}.faa -a pr_run/proteins/${i%%.fa}.faa;
		# run eggNOG
		cd /scratch/sogin/tmp.$JOB_ID/pr_run/proteins/
		python2 ~/tools/eggnog-mapper/emapper.py -i ${i%%.fa}.faa --output ${i%%.fa}_maNOG -m diamond --usemem --cpu 12;
done

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/edge_bins/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
