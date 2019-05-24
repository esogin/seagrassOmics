#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q@@himem
#
# Run Quast
#
#echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/metaquast/ /scratch/sogin/tmp.$JOB_ID
cd /scratch/sogin/tmp.$JOB_ID/
#
## Run QUAST 
#metaquast.py 3847_A.contigs.fasta  3847_C.contigs.fasta  3847_E.contigs.fasta 3847_F.contigs.fasta 3847_H.contigs.fasta coassembly.contigs.fa 3847_B.contigs.fasta  3847_D.contigs.fasta  3847_G.contigs.fasta  3847_I.contigs.fasta -o combined_report;
metaquast.py *.contigs.fa* -o combined_report;
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/metaquast/; #change dir according to project
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
