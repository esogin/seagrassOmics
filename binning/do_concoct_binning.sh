#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# Binning script with Coconct 1.0.0
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/concoct/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID
#
conda activate concoct_env
cut_up_fasta.py contigs.fa  -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed bam/*sorted.bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -t 48 -b result/
merge_cutup_clustering.py result/clustering_gt1000.csv > result/clustering_merged.csv
mkdir result/fasta_bins
extract_fasta_bins.py contigs.fa result/clustering_merged.csv --output_path result/fasta_bins
#
conda deactivate concoct_env
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/concoct/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
