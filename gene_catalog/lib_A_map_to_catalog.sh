#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
#Template script to map libraries to gene catalog
#notes: the libraries being mapped have not been error corrected with spades and are NOT kmer normalized, hence the nn director name (for no normalization)
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/cazymes/maps_nn/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
echo "map to catalog";
#
#cd /scratch/sogin/tmp.$JOB_ID/maps/
lib=3847_A;
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/qc/"$lib"_q2_ktrimmed.fq.gz ./;
bbmap.sh in=${lib}_q2_ktrimmed.fq.gz ref=genes_cazy_annotation.fa covstats="$lib"_covstats out="$lib".bam scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
rm *fq.gz
rm genes_cazy_annotation.fa
#rm -r ref/
#
mkdir /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/cazymes/maps_nn/"$lib"/
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/cazymes/maps_nn/"$lib"/;
#rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
