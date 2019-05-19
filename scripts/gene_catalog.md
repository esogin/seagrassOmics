Created a gene catalog from both assembled reads 

All contig files were pasted together (cat) and renamed with unique fasta headers

Create a new gene catalog from all assembled sequences 


1. Dereplicate sequences, cluster sequences, call genes, map sequences back to catalog
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
#$ -l h="gc-node-21"
#
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
#
echo "cluster sequences with mmseqs"
#
#
mkdir /scratch/sogin/tmp.$JOB_ID/clusters/
#vsearch --derep_fulllength allcontigs_renamed.fa --output clusters/rep_set.fa  --relabel group --relabel_keep
cd /scratch/sogin/tmp.$JOB_ID/clusters/
#mmseqs easy-linclust clusters/rep_set.fa clusters tmp    
#
#
echo "call genes"
#
#
cd /scratch/sogin/tmp.$JOB_ID/
mkdir /scratch/sogin/tmp.$JOB_ID/genes
prodigal -i clusters/clusters_rep_seq.fasta -o genes/coords.gbk  -a genes/orfs.faa  -p meta -d genes/genes.fa 
#
#
echo "map sequences to clustered genes" 
#
#
mkdir /scratch/sogin/tmp.$JOB_ID/maps
libs=$(echo 3847_{A..I})
for lib in $libs; do
	in1=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
	in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;;
bbmap.sh in=$in1 in2=$in2 ref=genes/genes.fa covstats=maps/"$lib"_covstats out=maps/"$lib"_mapped.bam scafstats=maps/"$lib"_scafstats.txt statsfile=maps/"$lib"_stderr;
done
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

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
gene_catalog=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/genes/gene_catalog.dmnd;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/libraries_3847/gene_catalog2/maps/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
mkdir /scratch/sogin/tmp.$JOB_ID/diamond/;
cd /scratch/sogin/tmp.$JOB_ID/diamond/;
cp $in ./ 
## DIAMOND BLASTX Search Against custom Protein DB
diamond blastx -d ${gene_catalog} -q ${lib}_q2_ktrimmed.fq.gz -o ${lib}_gene_catalog.m8 -p 24 -k 1 --id 50 -e 0.00001;
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/${lib}/gene_catalog2/maps/; 
rm /scratch/sogin/tmp.$JOB_ID/ -R;
echo "job finished: "
date
```



2. Annotate Sequeces with CAZYmes and Prokka 

```Bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate with CAZYs
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/annotate_me/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
#Annotate with HMM scan & cazy database 
python ~/tools/run_dbcan/run_dbcan.py orfs.faa protein --out_dir output_orfs --db_dir ~/tools/run_dbcan/db/
#
# Annotate with Prokka
 prokka --outdir prokka_result  --evalue 0.001 --metagenome  

#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/annotate_me/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

3. Get maps out and gene counts 

```Bash

```






Steps: 
1. call genes using prodigal
2. dereplciate and cluster sequences
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 64
#$ -V
#$ -q main.q@@himem
#
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"



mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/mmseqs/ /scratch/sogin/tmp.$JOB_ID
cd /scratch/sogin/tmp.$JOB_ID/

# call open reading frames 
#mkdir /scratch/sogin/tmp.$JOB_ID/orfs
#cd /scratch/sogin/tmp.$JOB_ID/orfs
#prodigal -i allcontigs_renamed.fa -o coords.gbk  -a orfs.faa  -p meta -d genes.fa;
# dereplicate ORFs
#vsearch --derep_fulllength genes.fa --output rep_set.fa  --relabel_ge1 --relabel_keep
#vsearch --threads 64 --log CLUSTER_CONTIGS --cluster_fast rep_set.faa --centroids centriods_orfs.faa --id 0.95 --consout concensus_orfs.faa --sizeout;
# Alternative clustering workflow using mmseqs
vsearch --derep_fulllength genes.fa --output rep_set.fa  --relabel group --relabel_keep
mmseqs easy-linclust rep_set.fa clusters tmp     
#
#
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/mmseqs/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

3. annotate sequences 
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
## Derepcliate & Annotate sequences
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/ /scratch/sogin/tmp.$JOB_ID
cd /scratch/sogin/tmp.$JOB_ID/
#
#vsearch --derep_fulllength allcontigs_renamed.fa --output rep_set.fa  --relabel group
#
# for some readson mmseqs only works on node 21
#
##mmseqs easy-linclust rep_set.fa clusters tmp   
prokka clusters_rep_seq.fasta --metagenome --outdir annotations
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```


Repeat annotation but with cazymes database and hmm scaner
```Bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q@@himem
#$ -l h="gc-node-21"
## Derepcliate & Annotate sequences
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/cazy_catalog/ /scratch/sogin/tmp.$JOB_ID
cd /scratch/sogin/tmp.$JOB_ID/
#
#Annotate with HMM scan & cazy database 
vsearch --derep_fulllength allcontigs_renamed.fa --output rep_set.fa  --relabel group --relabel_keep
mmseqs easy-linclust rep_set.fa clusters tmp     
prodigal -i clusters_rep_seq.fasta -o coords.gbk  -a orfs.faa  -p meta -d genes.fa;
hmmscan --domtblout proteins_annotated_cazy.out.dm cazy/dbCAN-fam-HMMs.txt orfs.faa > proteins_annotated_cazy.out;
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/cazy_catalog/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```


4. Create mapping script for cleaned libraries to be mapped back to clustered gene catalog

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
# Assembly script - maxbin
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/mapping/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/
#
libs=$(echo 3847_{A..I})
for lib in $libs; do
in1=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
in3=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers__unpaired.00.0_0.cor.fastq.gz;
  bbwrap.sh in=$in1,$in2,$in3 ref=genes.fa covstats=mapped/"$each"_covstats out=mapped/"$each"_mapped.bam scafstats=mapped/"$each"_scafstats.txt statsfile="$each"_stderr;
done
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/mapping/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

















