#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
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
lib=3847_I;
ref=/opt/extern/bremen/symbiosis/tools_HGV/adapters.fa;
coassembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/megahit/coassembly.contigs.fixed.fa;
assembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/assembly/spades/contigs.fasta
in=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/stats/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/

#fix fasta headers for down stream analysis
reformat.sh in=$assembly out=contigs.fixed.fasta addunderscore
# Map reads back to assembly
#
bbmap.sh ref=contigs.fixed.fasta in=$in in2=$in2 minid=0.98 fast=t covstats=${lib}_COV_STATS statsfile=${lib}_STATS outm=${lib}_mapped.bam
samtools sort -o ${lib}_mapped_sorted.bam ${lib}_mapped.bam
samtools index ${lib}_mapped_sorted.bam
#
# Map reads back to coassembly
bbmap.sh bbmap.sh ref=$coassembly in=$in in2=$in2 minid=0.98 fast=t covstats=${lib}_COV_STATS_coassembly statsfile=${lib}_STATS_coassembly outm=${lib}_mapped_coassembly.bam
samtools sort -o ${lib}_mapped_coassembly_sorted.bam ${lib}_mapped_coassembly.bam
samtools index ${lib}_mapped_coassembly_sorted.bam
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/stats;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
