#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
## Profile with anvio

echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/anvio/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
conda activate anvio5
#Reformat fasta file of contigs and remove contigs with length less then 1000
#Keep bin names the same so that we can import other mapping data onto contigs database at end
#anvi-script-reformat-fasta contigs.fa -o contigs-fixed.fa -l 1000 
#mv contigs-fixed.fa contigs.fa
#anvi-gen-contigs-database -f contigs.fa -o contigs.db -n "Sediment MG Database"
#anvi-run-hmms -c contigs.db --num-threads 48
#anvi-display-contigs-stats contigs.db  # only run if have ssh tunnel in place 
#anvi-setup-ncbi-cogs #do only once
#Annotate assembly with cog categories, may be useful later on
#anvi-run-ncbi-cogs -c contigs.db --num-threads 48
#Initialize BAM files
#bams=$(echo 3847_{A..I}_sorted.bam)
#for sample in $bams; do anvi-init-bam bam/$sample -o $sample; done
#Profile BAM files
bams=$(echo 3847_{A..I})
for i in $bams;do
	anvi-profile -i data/${i}_sorted.bam -c contigs.db --output-dir ${i}_profile --sample-name profile_"$i" -T 4;
done

conda deactivate
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/anvio/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
