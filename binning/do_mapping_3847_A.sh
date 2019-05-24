#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# Mapping script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
#lib=LIB;
lib=3847_A;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
cd /scratch/sogin/tmp.$JOB_ID/
#
# Map reads back to coassembly
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz ./
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz ./
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/contigs.fa ./
#
in1=${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
in2=${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
#
mkdir /scratch/sogin/tmp.$JOB_ID/bam/
bbmap.sh in=$in1 in2=$in2 ref=contigs.fa covstats="$lib"_COVSTATS out=bam/"$lib".bam scafstats="$each"_scafstats.txt statsfile="$each"_stderr;
samtools sort -o bam/${lib}_sorted.bam bam/"$lib".bam
samtools index bam/${lib}_sorted.bam 
#
#remove extra files
rm contigs.fa ${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz ${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maps/${lib}/;
#rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
echo "clean up scratch: sogin/tmp.$JOB_ID"
date
