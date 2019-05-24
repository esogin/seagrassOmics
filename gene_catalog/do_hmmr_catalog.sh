#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate with CAZYs
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/hmmr_annotation/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
#hmmr way
hmmscan -o cazy.orfs.out --tblout cazy.orfs.tblout --pfamtblout cazy.orfs.pfamblout --cpu 24 --domtblout cazy.orfs.out.dm caZy/dbCAN-fam-HMMs.txt orfs.faa;
hmmscan -o pfam.orfs.out --tblout pfam.orfs.tblout --pfamtblout pfam.orfs.pfamblout --cpu 24 --domtblout pfam.orfs.out.dm pfam/Pfam-A.hmm orfs.faa;
#./hmmscan-parser.sh cazy.orfs.out.dm > cazy_annotated_parsed.txt
#./hmmscan-parser.sh pfam.orfs.out.dm > pfam_annotated_parsed.txt
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/hmmr_annotation/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
