# Anotate select MAGs
## EM Sogin
## Created: May 22, 2019
## Updated: May 23, 2019


### 1. Annotate mags with antismash

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
## Annotate bins with prokka
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/bins_of_interest/bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
conda activate antismash
cd /scratch/sogin/tmp.$JOB_ID/bac/
for i in *fa; do 
	mkdir ${i%%.fa}
	prodigal -i $i -o ${i%%.fa}/${i%%.fa}_coords.gbk  -a ${i%%.fa}/${i%%.fa}_orfs.faa -d ${i%%.fa}/${i%%.fa}_genes.fa
 	antismash ${i%%.fa}/${i%%.fa}_orfs.fa -c 48 --taxon bacteria --outputfolder ${i%%.fa} --full-hmmer
done
rm *fa
#
conda deactivate
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/bins_of_interest/antismash/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Notes on AntiSMASH Analysis
1. The above script works but only at the moment with antismash 4.1, instead upload mags to webserver interface and save results for interperation later. 
(OK- only 25 mags, if more then consider reimplementing)

2. Initial impressions 
	2a. MAG metabat.232 contains 67 gene clusters that have antismash ID's domains, this is pretty high given that most MAGS (in general) have much lower or no antibiotic gene clusters. This mag is a Cellvibrionaceae and classified iwthin the genus teredinibacter, would be interesting to see how closely related it with other members of the group. However, elevated abundance really only higher inside one library within the meadow. Check 16S reads. 

	2b. In general, our mags have biosynthetic gene clusters for antiSMASH natural products. A lot are caracterized in the terpene and NRPS (non ribosomal peptide synthetase cluster) categories. Terpenes are interesting as they are largely produced by plants. **Might be worth while checking OTHER outside MAGS for gimilar gene clusters**

	2c. Idea: Check seagrass genome paper for terpenes in genome reptor, do we know if seagrasses are producing them? Has this been characterized. 
		- in genome paper suggests that terpenoid genes are reduced to only 2 genes, anything else making terpenoids? 


### 2. Annotate genomes in Patric

1. Sediment bins that are more abundant inside the seagrass meadows then outside were uploaded to patric on May 22, 2019. 

2. In order to annotate each bin in patric, we need to privide some level of classification. Try to use the GTBTK results to fill in this field and proceed to annotate all bins. Bin anntoation results stored in subfolder: annotation. 



### 3. Annotate genomes with prokka
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate bins with prokka
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/bins_of_interest/bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
cd /scratch/sogin/tmp.$JOB_ID/arc/
for i in *fa; do 
	prokka  --cpus 24 --kingdom Archaea --outdir "$i"_out $i; 
done
rm *fa
#
cd /scratch/sogin/tmp.$JOB_ID/bac/
for i in *fa; do 
	prokka --cpus 24 --outdir prokka_"$i"_result $i; 
done
rm *fa
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/bins_of_interest/prokka/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```


### 4. Annotate bins with COG categories 

steps: 
1. Call genes (prodigal)
2. use the eggNOG mapper to call anntoations and get OGs

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 12
#$ -V
#$ -q main.q
## Annotate bins with eggNOG
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/meadow_bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/bins/
bins=$(echo ls *fa)

cd /scratch/sogin/tmp.$JOB_ID/

for i in $bins; 
	do
		# run prodigial 
		prodigal -i bins/$i -o pr_run/genes/${i%%.fa}.faa -a pr_run/proteins/${i%%.fa}.faa;
		# run eggNOG
		cd /scratch/sogin/tmp.$JOB_ID/pr_run/proteins/
		python2 ~/tools/eggnog-mapper/emapper.py -i ${i%%.fa}.faa --output ${i%%.fa}_maNOG -m diamond --usemem --cpu 1;
done

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/meadow_bins/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```
Get cog results into useable format 
```bash
files=$(echo *emapper.annotations)
for i in $files; do cut -f 1,21,22 $i > ${i%%_maNOG.emapper.annotations}_COG.txt; done
#replace tabs with ; 

files=$(echo *COG.txt)
for i in $files; do  sed "s/\t/;/g" $i > ${i%%_COG.txt}_COG_fixed.txt; done

```

```R
files<-list.files(pattern="COG_fixed.txt", recursive=T)
# sep in to in-edge-out bins
in_files<-files[grep('meadow',files)]
df.in<-data.frame()
for (i in 1:length(in_files)){
df<-read.table(in_files[i],header=F, sep=';', quote="",fill=T)
df$binID<-in_files[i]
df$site<-'meadow'
df.in<-rbind(df.in, data.frame(df))
}

edge_files<-files[grep('edge',files)]
df.edge<-data.frame()
for (i in 1:length(edge_files)){
df<-read.table(edge_files[i],header=F, sep=';', quote="",fill=T)
df$binID<-edge_files[i]
df$site<-'edge'
df.edge<-rbind(df.edge, data.frame(df))
}

out_files<-files[grep('out',files)]
df.out<-data.frame()
for (i in 1:length(out_files)){
df<-read.table(out_files[i],header=F, sep=';', quote="",fill=T)
df$binID<-out_files[i]
df$site<-'out'
df.out<-rbind(df.out, data.frame(df))
}

cogs<-rbind(df.in, df.edge, df.out)
write.csv(cogs, file='bin_cog_annotation.csv')

```

### 5. Annotate with DBCan 
Lessons: need to have dbcan in directory or really have the file of interest in the DB can directory :eyeroll: 

All proteins were first called with prodigial and then run through DB can pipeline

```Bash 
#!/bin/bash
#
#$ -cwd 
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate bins run dbcan 
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 

rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/target_bins/out_bins/run_dbcan/ /scratch/sogin/tmp.$JOB_ID/
bins=$(echo *.faa)
for i in $bins; do 
 python run_dbcan.py $i protein --out_dir ${i%%.faa}_output --dia_cpu 24 --hmm_cpu 24 --hotpep_cpu 24;

 python run_dbcan.py 102_sub.contigs_proteins.faa protein --out_dir output_test --dia_cpu 24 --hmm_cpu 24 --hotpep_cpu 24 
done
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/target_bins/out_bins/run_dbcan/ 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Run DB can on all bins reguadless of quality/ habitat binning etc. Sort these bins out later
```Bash 
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q
## Annotate bins run dbcan 
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/ /scratch/sogin/tmp.$JOB_ID/

bins=$(echo ls *fa)
for i in $bins; 
	do
		# run prodigial 
		prodigal -i bins/$i -o genes/${i%%.fa}_genes.fa -a proteins/${i%%.fa}_proteins.faa;
done

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/ 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```


get data out of dbcan result
move results folder to new folder 
```r
files<-list.files(pattern='overview.txt', recursive=T)
dbcan_results<-data.frame()
for (i in 1:length(files)){
df<-read.table(files[i],sep='\t',header=T)
colnames(df)<-c('Gene.ID', 'HMMER', 'Hotpep' ,'DIAMOND', 'Signalp', 'num_tools')
df$bin<-files[i]
dbcan_results<-rbind(dbcan_results, df)
}
```
