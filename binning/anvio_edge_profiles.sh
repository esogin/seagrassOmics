#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
## Profile with anvio

echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/anvio/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID
conda activate anvio5
#Library A
anvi-profile -i data/3847_D_sorted.bam -c contigs.db --output-dir 3847_D_profile --sample-name profile_3847_D -T 1 --min-contig-length 5000;
anvi-profile -i data/3847_E_sorted.bam -c contigs.db --output-dir 3847_E_profile --sample-name profile_3847_E -T 1 --min-contig-length 5000;
anvi-profile -i data/3847_F_sorted.bam -c contigs.db --output-dir 3847_F_profile --sample-name profile_3847_F -T 1 --min-contig-length 5000;
conda deactivate
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/anvio/edge_profiles/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
