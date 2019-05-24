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
lib=3847_F;
#ref=/opt/extern/bremen/symbiosis/tools_HGV/adapters.fa;
#coassembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/megahit/coassembly.contigs.fixed.fa;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p;
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/3847_F/assembly/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID -p;
#
# Quality control filtering, trimming & normalize
#mkdir /scratch/sogin/tmp.$JOB_ID/qc/
#cd /scratch/sogin/tmp.$JOB_ID/qc/
#
#bbduk.sh ref=${ref} ktrim=l mink=11 hdist=1 in=/scratch/sogin/tmp.$JOB_ID/${lib}_R1.fastq.gz in2=/scratch/sogin/tmp.$JOB_ID/${lib}_R2.fastq.gz out=${lib}_ktriml.fq.gz;
#
#bbduk.sh ref=${ref} ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 in=${lib}_ktriml.fq.gz out=${lib}_q2_ktrimmed.fq.gz;
#
#bbnorm.sh -Xmx400g  in=${lib}_q2_ktrimmed.fq.gz out=${lib}_highfreq_kmers.fq.gz target=100 min=2;
#
# Run metaspades
#mkdir /scratch/sogin/tmp.$JOB_ID/assembly/
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/qc/${lib}_highfreq_kmers.fq.gz  /scratch/sogin/tmp.$JOB_ID -p;
spades.py -o spades_v --12 ${lib}_highfreq_kmers.fq.gz -t 48 -m 850 --phred-offset 33;
#
#fix fasta headers for down stream analysis
#reformat.sh in=contigs.fasta out=contigs.fixed.fasta addunderscore

#
# Map reads back to assembly
#mkdir /scratch/sogin/tmp.$JOB_ID/stats/
#cd /scratch/sogin/tmp.$JOB_ID/stats/
#bbmap.sh ref=../assembly/spades/contigs.fixed.fasta in=../spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz in2=../spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz minid=0.98 fast=t covstats=${lib}_COV_STATS statsfile=${lib}_STATS outm=${lib}_mapped.bam
#
# Map reads back to coassembly
#bbmap.sh ref=${coassembly} in=../spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz in2=../spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz minid=0.98 fast=t covstats=${lib}_COV_STATS statsfile=${lib}_STATS outm=${lib}_mapped_coassembly.bam
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
