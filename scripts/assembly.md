### Sediment metagenome assemblies
# EM Sogin
# Created Jan 14, 2019
# Updated March 28 2019

Sediment assembly script

Preform error correction on 2 TB server and run megahit to coassembly all libraries
```bash

# Quality control filtering and trimming
bbduk.sh ref=adapters.fa ktrim=l mink=11 hdist=1 in=3847_R1.fastq.gz in2=3847_R2.fastq.gz out=3847_ktriml.fq.gz;
bbduk.sh ref=adapters.fa ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 in=3847_ktriml.fq.gz out=3847_q2_ktrimmed.fq.gz;

# Normalize kmer spectra
bbnorm.sh -Xmx400g  in=3847_q2_ktrimmed.fq.gz out=3847_highfreq_kmers.fq.gz target=100 min=2;

# error correct with spades 
spades.py -o spades --12 3847_highfreq_kmers.fq.gz -t 64 -m 1.8T --phred-offset 33 --only-error-correction;

#Run megahit
megahit -1 spades/corrected/3847_highfreq_kmers_1.00.0_0.cor.fastq.gz -2 spades/corrected/3847_highfreq_kmers_2.00.0_0.cor.fastq.gz -r spades/corrected/3847_highfreq_kmers__unpaired.00.0_0.cor.fastq.gz -t 48 -o megahit --out-pre 3847 --k-min 21 --k-max 151 --k-step 10

# Map back to assembly
bbmap.sh ref=megahit/coassembly.contigs.fa in=corrected/3847_highfreq_kmers_1.00.0_0.cor.fastq.gz in2=corrected/3847_highfreq_kmers_2.00.0_0.cor.fastq.gz minid=0.98 fast=t covstats=3847_COV_STATS statsfile=3847_STATS outm=3847_mapped.sam
```

Template script for each library assembly 
```bash
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
lib=LIB;
ref=/opt/extern/bremen/symbiosis/tools_HGV/adapters.fa;
coassembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/library_3847/coassembly/megahit/coassembly.contigs.fixed.fa;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/ /scratch/sogin/tmp.$JOB_ID; 
#
# Quality control filtering, trimming & normalize
mkdir /scratch/sogin/tmp.$JOB_ID/qc/
cd /scratch/sogin/tmp.$JOB_ID/qc/
#
bbduk.sh ref=${ref} ktrim=l mink=11 hdist=1 in=/scratch/sogin/tmp.$JOB_ID/${lib}_R1.fastq.gz in2=/scratch/sogin/tmp.$JOB_ID/${lib}_R2.fastq.gz out=${lib}_ktriml.fq.gz;
#
bbduk.sh ref=${ref} ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 in=${lib}_ktriml.fq.gz out=${lib}_q2_ktrimmed.fq.gz;
#
bbnorm.sh -Xmx400g  in=${lib}_q2_ktrimmed.fq.gz out=${lib}_highfreq_kmers.fq.gz target=100 min=2;
#
# Run metaspades
mkdir /scratch/sogin/tmp.$JOB_ID/assembly/
cd /scratch/sogin/tmp.$JOB_ID/assembly/
spades.py -o spades --12 ../qc/${lib}_highfreq_kmers.fq.gz -t 24 -m 850 --phred-offset 33;
#
#fix fasta headers for down stream analysis
reformat.sh in=contigs.fasta out=contigs.fixed.fasta addunderscore
#
#
# Map reads back to assembly
mkdir /scratch/sogin/tmp.$JOB_ID/stats/
cd /scratch/sogin/tmp.$JOB_ID/stats/
bbmap.sh ref=../assembly/spades/contigs.fixed.fasta in=../spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz in2=../spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz minid=0.98 fast=t covstats=${lib}_COV_STATS statsfile=${lib}_STATS outm=${lib}_mapped.bam
#
# Map reads back to coassembly
bbmap.sh ref=${coassembly} in=../spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz in2=../spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz minid=0.98 fast=t covstats=${lib}_COV_STATS statsfile=${lib}_STATS outm=${lib}_mapped_coassembly.bam
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
```
Mapping above failed, repeat with code below and new template script: mapping

```bash
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
lib=LIB;
ref=/opt/extern/bremen/symbiosis/tools_HGV/adapters.fa;
coassembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/megahit/coassembly.contigs.fixed.fa;
assembly=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/contigs.fasta
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
bbmap.sh ref=$coassembly in=$in in2=$in2 minid=0.98 fast=t covstats=${lib}_COV_STATS_coassembly statsfile=${lib}_STATS_coassembly outm=${lib}_mapped_coassembly.bam

bbmap.sh ref=$ref in=$in1 in2=$in2 minid=0.98 fast=t covstats=COV_STATS statsfile=STATS outm=mapped.bam

samtools sort -o ${lib}_mapped_coassembly_sorted.bam ${lib}_mapped_coassembly.bam
samtools index ${lib}_mapped_coassembly_sorted.bam
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/stats;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
```

run metaquost
```bash
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
metaquast.py 3847_A.contigs.fasta  3847_C.contigs.fasta  3847_E.contigs.fasta 3847_F.contigs.fasta 3847_H.contigs.fasta coassembly.contigs.fa 3847_B.contigs.fasta  3847_D.contigs.fasta  3847_G.contigs.fasta  3847_I.contigs.fasta -o combined_report;
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/metaquast/; #change dir according to project
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
```






