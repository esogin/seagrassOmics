#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/cazymes/maps/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
echo "map to catalog";
#
#cd /scratch/sogin/tmp.$JOB_ID/maps/
libs=$(echo 3847_{A..I})
for lib in $libs; do
	cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz ./;
	cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz ./;
	bbmap.sh in=${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz in2=${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz ref=genes_cazy_annotation.fa covstats="$lib"_covstats out="$lib".bam scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
done
rm *fastq.gz
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/cazymes/maps/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
