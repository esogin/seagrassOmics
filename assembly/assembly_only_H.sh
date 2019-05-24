#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
#
# Assembly script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
lib=3847_H;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/ /scratch/sogin/tmp.$JOB_ID; 
#
# Quality control filtering, trimming & normalize
cd /scratch/sogin/tmp.$JOB_ID/
spades.py -o spades_v2 --12 ${lib}_highfreq_kmers.fq.gz -t 48 -m 850 --phred-offset 33;
#fix fasta headers for down stream analysis
reformat.sh in=spades_V2/contigs.fasta out=spades_v2/contigs.fixed.fasta addunderscore
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
