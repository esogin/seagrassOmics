#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
#$ -l h="gc-node-21"
#
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/prodigial/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
#
#echo "cluster sequences with mmseqs"
#
#
#mkdir /scratch/sogin/tmp.$JOB_ID/clusters/
#vsearch --derep_fulllength allcontigs_renamed.fa --output clusters/rep_set.fa  --relabel group --relabel_keep
#cd /scratch/sogin/tmp.$JOB_ID/clusters/
#mmseqs easy-linclust clusters/rep_set.fa clusters tmp
#
#
echo "call genes"
#
#
#cd /scratch/sogin/tmp.$JOB_ID/
mkdir /scratch/sogin/tmp.$JOB_ID/genes
prodigal -i clusters_rep_seq.fasta -o genes/coords.gff  -a genes/orfs.faa  -p meta -d genes/genes.fa -f gff; 
#
#
#echo "map sequences to clustered genes" 
#
#
#mkdir /scratch/sogin/tmp.$JOB_ID/maps
#libs=$(echo 3847_{A..I})
#for lib in $libs; do
#	in1=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
#	in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
#bbmap.sh in=$in1 in2=$in2 ref=genes.fa covstats="$lib"_covstats out="$lib".bam scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
#done
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/prodigial/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
