#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
#
# Assembly script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
lib=3847_B;
coassembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/megahit/coassembly.contigs.fixed.fa;
in=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/stats/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/

# Map reads back to coassembly
bbmap.sh ref=$coassembly in=$in in2=$in2 minid=0.98 fast=t covstats=${lib}_COV_STATS_CO statsfile=${lib}_STATS_CO outm=${lib}_mapped_coassembly.bam
samtools sort -o ${lib}_mapped_coassembly_sorted.bam ${lib}_mapped_coassembly.bam
samtools index ${lib}_mapped_coassembly_sorted.bam
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/stats;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
