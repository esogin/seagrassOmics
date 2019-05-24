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
lib=3847_H;
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
cp ${in} ./;
gunzip ${lib}_q2_ktrimmed.fq.gz;
kaiju -t ${kaiju_nodes} -f ${kaiju_fmi} -i ${lib}_q2_ktrimmed.fq -z 8 -a "greedy" -e 1 -o ${lib}_kaiju;
kaijuReport -t ${kaiju_nodes} -n ${kaiju_names} -i ${lib}_kaiju -r genus -p -m 0.1 -o ${lib}_kaiju.out.summary;
kaiju2krona -t ${kaiju_nodes} -n ${kaiju_names} -i ${lib}_kaiju -o ${lib}.krona;
perl /home/sogin/tools/KronaTools/scripts/ImportText.pl -o ${lib}_classification.html ${lib}.krona;
rm *.fq
rsync -a /scratch/sogin/tmp.$JOB_ID /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/read_based/taxonomy/; 
echo "job finished: "
date
