#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# binning script with maxbin2
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maxbin/ /scratch/sogin/tmp.$JOB_ID;
cd /scratch/sogin/tmp.$JOB_ID/
#
# Run maxbin2
run_MaxBin.pl -thread 48  -contig contigs.fa -out maxbin -abund 3847_A_COVSTATS_forMaxBin.txt -abund2 3847_B_COVSTATS_forMaxBin.txt -abund3 3847_C_COVSTATS_forMaxBin.txt -abund4 3847_D_COVSTATS_forMaxBin.txt -abund5 3847_E_COVSTATS_forMaxBin.txt -abund6 3847_F_COVSTATS_forMaxBin.txt -abund7 3847_G_COVSTATS_forMaxBin.txt -abund8 3847_H_COVSTATS_forMaxBin.txt -abund9 3847_I_COVSTATS_forMaxBin.txt
rsync -a /scratch/sogin/tmp.$JOB_ID /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maxbin/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
