# Sediment metagenome assemblies
## EM Sogin
## Created Jan 14, 2019
## Updated March 28 2019

Sediment assembly script
Description: this script seeks to assembly sediment metagenomic libraries using both servers in Cologne and in house at the MPI. The coassembly error correction needs to be done using at 2 TB memory server at the MPI. Megahit was used to assemnbly the libraries. For single library assemblies, the high mem servers in Cologne were sufficient and metaspades assembled the libraries in approximately 1 day/ library. 

## 1. Assembly libries
### Coassembly
Forward and reverse reads were concatinated across all 9 libraries. Preform error correction on 2 TB server and run megahit to coassembly all libraries.
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
### Library assembly
These are listed as bash qsub scripts for specific parameters in how they were run. 
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
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/ /scratch/sogin/tmp.$JOB_ID; 
#
# Quality control filtering, trimming & normalize with Kmer counts
#
mkdir /scratch/sogin/tmp.$JOB_ID/qc/
cd /scratch/sogin/tmp.$JOB_ID/qc/
bbduk.sh ref=${ref} ktrim=l mink=11 hdist=1 in=/scratch/sogin/tmp.$JOB_ID/${lib}_R1.fastq.gz in2=/scratch/sogin/tmp.$JOB_ID/${lib}_R2.fastq.gz out=${lib}_ktriml.fq.gz;
bbduk.sh ref=${ref} ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 in=${lib}_ktriml.fq.gz out=${lib}_q2_ktrimmed.fq.gz;
bbnorm.sh -Xmx400g  in=${lib}_q2_ktrimmed.fq.gz out=${lib}_highfreq_kmers.fq.gz target=100 min=2;
#
# Run metaspades
mkdir /scratch/sogin/tmp.$JOB_ID/assembly/
cd /scratch/sogin/tmp.$JOB_ID/assembly/
spades.py -o spades --12 ../qc/${lib}_highfreq_kmers.fq.gz -t 24 -m 850 --phred-offset 33;
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
```
##2. Check assembly stats
Map reads back to individual library assemblies 
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

#Fix fasta headers for down stream analysis of assemblies
reformat.sh in=$assembly out=contigs.fixed.fasta addunderscore
#
# Map reads back to assembly
#
bbmap.sh ref=contigs.fixed.fasta in=$in in2=$in2 minid=0.98 fast=t covstats=${lib}_COV_STATS statsfile=${lib}_STATS outm=${lib}_mapped.bam
samtools sort -o ${lib}_mapped_sorted.bam ${lib}_mapped.bam
samtools index ${lib}_mapped_sorted.bam
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/stats;
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
date
```
### Run metaquost for stats on assemblies
First move all assemblies to same folder to facilitate analysis
Results stored in html file in same folder
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

##END