#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
# Binning script with metabat
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/metabat/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/
# run metabat2
jgi_summarize_bam_contig_depths --minContigLength 1000 --outputDepth depth.txt --pairedContigs paired.txt bam/*sorted.bam 
metabat2 -i contigs.fa -a depth.txt -o result/metabat -v 
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/metabat/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
