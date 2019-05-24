#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q low.q
#
# Read based analysis - Function
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
lib=3847_B;
in=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/qc/${lib}_q2_ktrimmed.fq.gz;
sugars=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/sugar_proteins/sugars_uniport_prok.dmnd;
cherrys=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/cherrys/cherry_and_anticherry.dmnd;
cazy=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/cazy/caZy.dmnd;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/read_based/function/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
mkdir /scratch/sogin/tmp.$JOB_ID/diamond/;
cd /scratch/sogin/tmp.$JOB_ID/diamond/;
cp $in ./ 
#bbmerge.sh in=${lib}_q2_ktrimmed.fq.gz out=${lib}_merged.fq.gz ihist=${lib}_ihist.txt;
#reformat.sh in=${lib}_merged.fq.gz reads=10000000 out=${lib}_even_read_counts.fq.gz
#Run Humann2
#humann2 --input ${lib}_merged.fq.gz --output ${lib}_humann2
#Run sugar clasifications 
## DIAMOND BLASTX Search Against custom Protein DB
diamond blastx -d ${cherrys} -q ${lib}_q2_ktrimmed.fq.gz -o ${lib}_matches_cherry.m8 -p 24 -k 1 --id 50 -e 0.00001;
diamond blastx -d ${cazy} -q ${lib}_q2_ktrimmed.fq.gz -o ${lib}_matches_cazy.m8 -p 24 -k 1 --id 50 -e 0.00001;
diamond blastx -d ${sugars} -q ${lib}_q2_ktrimmed.fq.gz -o ${lib}_matches_sugars.m8 -p 24 -k 1 --id 50 -e 0.00001;

#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/read_based/function/; 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
