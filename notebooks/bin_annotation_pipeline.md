#bin annotation pipeline
## EM Sogin 
## Created July 17, 2019
## Updated July 17 2019

Redoing bin annotation pipeline to clarify and apply to all high QC bins, not just the target bins. 


base direcotry: 
/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/

result of dastool combination of bins

### Idnetification of habitat associated bins

1. Used visualization pipeline/ analysis notebook to determine which habitats the high quality bins were associated with. Some refinement and identification of bins that are inside, at the edge, outside the meadow. As well as identification of singleton bins. 

Habitat assocation is supported with a mixture of statistics, coverage values and visualization approaches. If the bin was only associated iwth a single sample, it was considered noise in the analysis and removed. If the average coverage across libraries was 1 it was considered present in the habitat. 


Results: based on presence (avg fold value > 1)
| habitat | number of bins|
|:-------:|:-------------:|
| All|14|
|In|16|
|Edge|21|
|Out|16|
|In_Edge|2|
|In_Out|5|
|Edge_Out|34|
|Single libraries|6|

Results: based on statistics (HSD test mean coverage across contigs ~ habitat) & visualization 

| habitat | number of bins|
|:-------:|:-------------:|
| All|11|
|In|17|
|Edge|24|
|Out|20|
|In_Edge|2|
|In_Out|3|
|Edge_Out|31|
|Single libraries|6|



**For comparison across libraries in annotation use the statistics plus visulization resutls in order to help determine which habitat bins belong to.**

### Gene calling

2. Use prodigal to call genes of high QC libraries

gene calling doing in prodigal direcorty

my qsub script

```bash
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
mkdir /scratch/sogin/tmp.$JOB_ID/genes/
mkdir /scratch/sogin/tmp.$JOB_ID/proteins/
cd /scratch/sogin/tmp.$JOB_ID/bins/ 
bins=$(echo *fa)
for i in $bins;
    do
# run prodigial 
       prodigal -i $i -o ../genes/${i%%.fa}_coords.gbk  -a ../proteins/${i%%.fa}_orfs.faa -d ../genes/${i%%.fa}_genes.fa;
done

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/ 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```  

### Annotation 

annotation done in highqc directory of prodigal directory

3. run dbcan 

copy all protein files to dbcan folder and run qsub script

```bash
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/run_dbcan/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID
bins=$(echo *.faa)
ls
echo $bins
for i in $bins; do
python run_dbcan.py $i protein --out_dir output_${i%%.faa} --dia_cpu 24 --hmm_cpu 24 --hotpep_cpu 24;
done
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/run_dbcan/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Use R to parse through dbcan results 
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
write.csv(dbcan_results, 'dbcan_results_July182019.csv')
```


4. run antismash 
For antismash, copy only bacteria genomes over into folder

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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/antismash/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
conda activate antismash
cd /scratch/sogin/tmp.$JOB_ID/bac/
for i in *fa; do 
    mkdir ${i%%.fa}
    prodigal -i $i -o ${i%%.fa}/${i%%.fa}_coords.gbk  -a ${i%%.fa}/${i%%.fa}_orfs.faa -d ${i%%.fa}/${i%%.fa}_genes.fa
    antismash ${i%%.fa}/${i%%.fa}_genes.fa -c 48 --taxon bacteria --input-type nucl --outputfolder ${i%%.fa} --full-hmmer --smcogs --transatpks_da
done
#
conda deactivate
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/antismash/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Get out antismash results 
```bash
find ./ -size  0 -print0 |xargs -0 rm --

files=$(ls)
for i in $files;do  cut -f 1,3,4,5 $i/geneclusters.txt > $i/geneclusters_fixed.txt;done
```

```r
files<-list.files(pattern='geneclusters_fixed.txt', recursive=T)
remove<-c(19,28,40,42,75,77,79,83,92,95,96,100)

asmsh<-data.frame()
for (i in 1:length(files)){
	tmp<-read.table(files[i],sep='\t')
	asmsh<-rbind(asmsh, data.frame(bin=files[i],tmp))
}

write.csv(asmsh, 'combined_antismash_annotations.csv')
```



5. run prokka 

copy bins over to prokka folder within the prodigial high qc bins directory
split bewteen archaea and bacteria 
run qsub script

qsub script
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/prokka/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
cd /scratch/sogin/tmp.$JOB_ID/arc/
for i in *fa; do
    prokka  --cpus 24 --kingdom Archaea --outdir "$i"_out $i;
done
#
cd /scratch/sogin/tmp.$JOB_ID/bac/
for i in *fa; do 
    prokka --cpus 24 --outdir prokka_"$i"_result $i;
done
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/prokka/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

6. run eggnog

Copy over proteins from prodigial run into eggnog folder 

run qsub script 

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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/eggnog/ /scratch/sogin/tmp.$JOB_ID/

cd /scratch/sogin/tmp.$JOB_ID/bins/
bins=$(echo *faa)

for i in $bins; do
        python2 ~/tools/eggnog-mapper/emapper.py -i $i --output ${i%%.fa}_maNOG -m diamond --usemem --cpu 12;
    done

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/prodigial/hiqc_bins/eggnog/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```



