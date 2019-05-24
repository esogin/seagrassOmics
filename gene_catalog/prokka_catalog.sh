#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
#
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/forProkka/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
prokka clusters_rep_seq.fasta --outdir result --metagenome --cpus 48 --mincontiglen 500
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/forProkka/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
