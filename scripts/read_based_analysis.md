#Read based analysis
##E. Sogin 
## March 29 2019

Script is set up for an alternative pipeline that takes reads and maps them to taxonomic classifiers and functional genes with diamond. Run script as taxonomy script and classification script for ease of use. 


```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
#
# Read based analysis - Taxonomy
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
lib=LIB;
in=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/qc/${lib}_q2_ktrimmed.fq.gz;

mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/read_based/taxonomy/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/

#Run phyloFlash
mkdir phyloFlash/;
cd phyloFlash/;
    phyloFlash.pl -read1 ${in} -Lib ${lib}_class -interleaved -CPUs 24 -zip -log; #class level taxonomy 
    phyloFlash.pl -read1 ${in} -Lib ${lib}_phyla -interleaved -taxlevel 3 -CPUs 24 -zip -log; 
    phyloFlash.pl -read1 ${in} -Lib ${lib}_order -interleaved -taxlevel 5 -CPUs 24 -zip -log; 
	phyloFlash.pl -read1 ${in} -Lib ${lib}_family -interleaved -taxlevel 6 -CPUs 24 -zip -log;

#Run Kaiju
mkdir /scratch/sogin/tmp.$JOB_ID/kaiju/
cd /scratch/sogin/tmp.$JOB_ID/kaiju/
kaiju_nodes=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/KAIJU/nodes.dmp;
kaiju_fmi=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/KAIJU/kaiju_db_nr_euk.fmi;
kaiju_nodes=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/KAIJU/nodes.dmp;
kaiju_names=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/databases/KAIJU/names.dmp;
#
gunzip ${lib}_q2_ktrimmed.fq.gz;
kaiju -t ${kaiju_nodes} -f ${kaiju_fmi} -i ${lib}_q2_ktrimmed.fq -z 8 -a "greedy" -e 1 -o ${lib}_kaiju;
kaijuReport -t ${kaiju_nodes} -n ${kaiju_names} -i ${lib}_kaiju -r genus -p -m 0.1 -o ${lib}_kaiju.out.summary;
kaiju2krona -t ${kaiju_nodes} -n ${kaiju_names} -i ${lib}_kaiju -o ${lib}.krona;
perl /home/sogin/tools/KronaTools/scripts/ImportText.pl -o ${lib}_classification.html ${lib}.krona;
rm *.fq
rsync -a /scratch/sogin/tmp.$JOB_ID /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/read_based/taxonomy/; 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Functional read classification
```bash
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
lib=LIB;
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
```