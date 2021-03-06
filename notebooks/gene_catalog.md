# Gene Catalog Approach 
## EM Sogin 
## Update: May 19, 2019

Created a gene catalog from assembled reads (both coassembly and individual assembly)

All contig files were pasted together (cat) and renamed with unique fasta headers
Create a new gene catalog from all assembled sequences 


## 1. Dereplicate sequences, cluster sequences, call genes, map sequences back to catalog

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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/maps/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
#
echo "cluster sequences with mmseqs"
#
#
mkdir /scratch/sogin/tmp.$JOB_ID/clusters/
vsearch --derep_fulllength allcontigs_renamed.fa --output clusters/rep_set.fa  --relabel group --relabel_keep
cd /scratch/sogin/tmp.$JOB_ID/clusters/
mmseqs easy-linclust clusters/rep_set.fa clusters tmp
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
    in2=/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
bbmap.sh in=$in1 in2=$in2 ref=genes.fa covstats="$lib"_covstats out="$lib".bam scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
done
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/maps/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```
## 2. Annotate Sequeces with CAZYmes (using hmmscan) and Prokka 
Prokka annotation will  likely take a long time to complete. Keep running in background.

Annotate using prokka
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q@@himem
#
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/forProkka/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
prokka clusters_rep_seq.fasta --outdir result --metagenome --cpus 48 --mincontiglen 500
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/forProkka/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Annotate using hmmer
```bash
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/cazymes/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
hmmscan --domtblout cazy.out.dm --cpus 24 db/dbCAN-fam-HMMs.txt orfs.faa > cazy.out 
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/cazymes/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```




Annotation using hmmer tools and run db can is more challening then initally considered. 

For a new appraoch, take gene catalog and first preform a translated blast search against cazymes database using diamond blastX. 
Contigs positive for cazyme hits, retain in dataset and use DB can (stand along on workstation) or on web interface to confirm diamond blastX hits. Diamond cut off use: e-value 10^-6 with 50% identity of AA sequences. 

```bash
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/diamond/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
diamond blastp --threads 24 -d ./db/CAZy.dmnd -q orfs.faa --evalue 0.000001 --id 50 --max-hsps 1 
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/diamond/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Check output -> evalues and % ID's correct? 
Take contigs with positive hits and subset orfs, run with DBcan


http://bcb.unl.edu/dbCAN2/download/



hmmr way --> old code
```bash
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/hmmr_annotation/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
#hmmr way
hmmscan -o cazy.orfs.out --tblout cazy.orfs.tblout --pfamtblout cazy.orfs.pfamblout --cpu 24 --domtblout cazy.orfs.out.dm caZy/dbCAN-fam-HMMs.txt orfs.faa;
hmmscan -o pfam.orfs.out --tblout pfam.orfs.tblout --pfamtblout pfam.orfs.pfamblout --cpu 24 --domtblout pfam.orfs.out.dm pfam/Pfam-A.hmm orfs.faa;
#./hmmscan-parser.sh cazy.orfs.out.dm > cazy_annotated_parsed.txt
#./hmmscan-parser.sh pfam.orfs.out.dm > pfam_annotated_parsed.txt
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/hmmr_annotation/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

This might be a good perl script to get hmmscan data out to a gff file
https://github.com/tomdeman-bio/HMMer-scripts-
https://mgkit.readthedocs.io/en/0.3.4/scripts/hmmer2gff.html





## 3. Generate featureCounts hit stats for BAM files



## 3. Use R helper script to combine results
prep scaff stats files for R import
```bash
cd /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/maps/
for i in $scafstats; do sed "s/ # /_/g" $i > fixed_"$i"; done
```
run commands in R
```r
## Combine read counts into single file for profiles
files<-list.files(pattern="fixed.*.txt")

df<-read.table(files[1],sep="\t", header=F)
colnames(df)<-c('name','%unambiguousReads',	'unambiguousMB','%ambiguousReads','ambiguousMB','unambiguousReads','ambiguousReads',	'assignedReads','assignedBases')
count<-data.frame(df$name, as.numeric(df$assignedReads))

colnames(count)<-c('name',paste0(files[1],sep=""))
#count$assignedReads<-as.vector(count$assignedReads)

for (i in 2:length(files)){
df2<-read.table(files[i],sep="\t", header=F)
colnames(df2)<-c('name','%unambiguousReads','unambiguousMB','%ambiguousReads','ambiguousMB','unambiguousReads','ambiguousReads',	'assignedReads','assignedBases')
count2<-data.frame(df2$name, as.numeric(df2$assignedReads))
colnames(count2)<-c('name',paste0(files[i],'_count',sep=""))
count<-merge(count, count2, by="name", all=T)
}

count<-count[is.na(count)]<-0
count$Total<-rowSums(count[,grep("fixed",colnames(count))])
range(count$Total)
summary(count$Total)

colSums(count[,grep("fixed",colnames(count))]) #between 49,820,708 and 101,056,230 reads maped to catalog

count.abund<-count[count$Total > 10, ]


## Add in gene identifiers to profile table
annotations<-read.table('~/extern/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/annotate_me/output_orfs/h.out_modified.csv',sep=',',header=T)
colnames(annotations)[1]<-'domain'

head(count.abund)
colnames(count.abund)[1]<-'name'
count.abund$name<-as.vector(count.abund$name)
count.abund$name2<-paste(sapply(strsplit(count.abund$name,"_"),'[',1),sapply(strsplit(count.abund$name,"_"),'[',2),sep="_")
count.abund$DomainID<-annotations[match(count.abund$name2,annotations$query_name),'domain']

# write to output
write.csv(count.abund,'gene_counts_cazyAnnot_May52019.csv')
```

## end















