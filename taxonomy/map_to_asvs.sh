!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 12
#$ -V
#$ -q main.q
#
# Mapping script - clean libraries to ASV file of 16S tags
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/taxonomy/map_to_16S_pb/maps/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/
#
ref=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/taxonomy/map_to_16S_pb/meadow_asvs.fa;
LIBS=$(echo 3847_{A..I})
for lib in $LIBS; do 
	in1=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
	in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
	bbmap.sh in=$in1 in2=$in2 ref=$ref minid=0.97 covstats="$lib"_COVSTATS scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
done
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/taxonomy/map_to_16S_pb/maps/
rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
echo "clean up scratch: sogin/tmp.$JOB_ID"
date
