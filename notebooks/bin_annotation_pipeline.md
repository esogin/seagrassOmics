#bin annotation pipeline
## EM Sogin 
## Created July 17, 2019
## Updated August 19, 2019

Redoing bin annotation pipeline to clarify and apply to all high QC bins, not just the target bins. 


base direcotry: 
/opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool_finalized_bins/result/

result of dastool combination of bins

Ran trna search across bins to help futher identify bins that are more complete. Went from 114 bins to 126 bins. Not a huge increase but still nice to have some more bins (~ 10% more)
Because of this, I am re-running the annotation and habtiat classificaiton pipelines previously estabolished. 

### Gene calling & annotation of high qc bins

1. Use prodigal to call genes of high QC libraries

my qsub script
Really important: the run_dbcan script must be in same folder as the proteins, see qsub script

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q

echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 

rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/final/high_qc_bins/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/

eggnog_path=/opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper
path_to_r_parsing_script=/opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper/parse_eggnog.R

# Run prodigal across all bins
echo "running prodigal"
mkdir /scratch/sogin/tmp.$JOB_ID/prodigal/
mkdir /scratch/sogin/tmp.$JOB_ID/prodigal/genes
mkdir /scratch/sogin/tmp.$JOB_ID/prodigal/coords
mkdir /scratch/sogin/tmp.$JOB_ID/prodigal/proteins

cd /scratch/sogin/tmp.$JOB_ID/bins/ 
bins=$(echo *fa)

for i in $bins;
    do
# run prodigial 
       prodigal -i $i -o ../prodigal/coords/${i%%.fa}_coords.gbk  -a ../prodigal/proteins/${i%%.fa}_orfs.faa -d ../prodigal/genes/${i%%.fa}_genes.fa;
done

# Run eggnog 
echo "running eggnog"

mkdir /scratch/sogin/tmp.$JOB_ID/eggnog/
cd /scratch/sogin/tmp.$JOB_ID/prodigal/proteins
bins=$(echo *faa)

for i in $bins; do
        python2 $eggnog_path/emapper.py -i $i --output ../../eggnog/${i%%.fa}_maNOG -m diamond --usemem --cpu 12;
    done

echo "running parsing script"
# Select emapper columns of interest 
# 1.  query_name, 2. seed eggNOG ortholog, 3. seed ortholog evalue, 7. Gene Ontology terms
# 8. EC number, 9. KEGG_Ko, 10. KEGG_Pathway, 11. KEGG_Module, 16. CAZy, 21. COG Functional Category
# see https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2 for other columns that you could use

mkdir /scratch/sogin/tmp.$JOB_ID/eggnog/new_result/
cd /scratch/sogin/tmp.$JOB_ID/eggnog/
emapper_results=$(ls *emapper.annotations)
for i in $emapper_results; do cut -f 1,2,3,7,8,9,10,11,16,21 $i > new_result/${i}_simplified; done

cd /scratch/sogin/tmp.$JOB_ID/eggnog/new_result/
Rscript $path_to_r_parsing_script 

echo "running run_db_can"
cd /scratch/sogin/tmp.$JOB_ID/run_dbcan
mkdir /scratch/sogin/tmp.$JOB_ID/run_dbcan/result
bins=$(echo proteins/*.faa)

for i in $bins; do
python2 run_dbcan.py $i protein --out_dir output_${i%%.faa} --dia_cpu 24 --hmm_cpu 24 --hotpep_cpu 24;
done

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/final/high_qc_bins/

rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```  

### ID of habitat associated bins

1. Use gb-tools to generate coverage statistics and create plots of individual bins 
bin list is a list of bins with bin file name and bin name

```bash
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i bin.list -o all_bins.tsv
```

```r
## Generate combined covstats file for analysis 
files<-list.files(pattern='forR',recursive=T, full.names=T)
LIBS<-c('3847_A','3847_B','3847_C','3847_D','3847_E','3847_F','3847_G','3847_H','3847_I')

#Read in covstats and paste into dataframe
covstats<-data.frame()
for (i in 1:length(files)){
lib_import<-read.table(files[i],header=T)
lib_import<-lib_import[which(lib_import$Length > 1000),]
lib_import$library<-LIBS[i]
covstats<-rbind(covstats, data.frame(lib_import))
}
# read bin files of meadow bins 
bins<-read.table('all_bins.tsv')
colnames(bins)<-c('name','contigID')
#match to names
covstats$name<-bins[match(covstats$ID,bins$contigID),'name'] 
covstats2<-covstats[complete.cases(covstats),]

#Save for plotting
write.csv(covstats2,'August2019_annotated_covstats_all_bins.csv')

### visualization with gbtools
library(gbtools)
d <- gbt(covstats=c("covstats/3847_A_COVSTATS", "covstats/3847_B_COVSTATS", "covstats/3847_C_COVSTATS","covstats/3847_D_COVSTATS","covstats/3847_E_COVSTATS", "covstats/3847_F_COVSTATS","covstats/3847_G_COVSTATS", "covstats/3847_H_COVSTATS","covstats/3847_I_COVSTATS"))
summary(d)

# Create plots for all lists --> inspect plots
pdf(file='gc_coverage_plots.pdf',width=15, height=15)
par(mfrow=c(3,3))
    for(i in 1:9){
        multiBinPlot(d, bins=d.bins, binNames=names(d.bins),slice=i,cutoff=10000, main=LIBS[i],legend=T)
    }
dev.off()


```

2. import data to Rnotebook for visualization and statistical anlayses 
R notebook is called: August 19 habitat bin assignmen.html

Bins were determined to be sig different between habitats based on an aov model and mean coverage of contigs > 5000 bp in length.
If model did not pick up sig results, the gc-coverage plots were checked to determine which habitat the bin exclusively occured in. 
This was because of our low sample size.

results

| habitat | number of bins|
|:-------:|:-------------:|
|Cosmopolatin|24|
|In|15|
|Edge|20|
|Out|21|
|in-edge|2|
|in-out|1|
|non-meadow|35|
|Single libraries|7|

Good news: similar number of bins bewteen habitats, 15-20-21 (in - edge - out)

Include non-meadow and cosmopolatin bins in futher analyses as well. 




---------------------------------------------------------------
---------------------------------------------------------------
## OLD CODE & Analyses 
---------------------------------------------------------------
---------------------------------------------------------------


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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/final/hiqc_bins_test/ /scratch/sogin/tmp.$JOB_ID
eggnog_path=/opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper

#R script will need to be modified if you want different columns from emapper in the cut command below
path_to_r_parsing_script=/opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper/parse_eggnog.R

mkdir /scratch/sogin/tmp.$JOB_ID/eggnog/
cd /scratch/sogin/tmp.$JOB_ID/proteins/
bins=$(echo *faa)

for i in $bins; do
        python2 $eggnog_path/emapper.py -i $i --output ../eggnog/${i%%.fa}_maNOG -m diamond --usemem --cpu 12;
    done

## Note: either run following commands in current datafolder or combine results from eggnog mapper with previous eggnog results.

# Select emapper columns of interest 
# 1.  query_name, 2. seed eggNOG ortholog, 3. seed ortholog evalue, 7. Gene Ontology terms
# 8. EC number, 9. KEGG_Ko, 10. KEGG_Pathway, 11. KEGG_Module, 16. CAZy, 21. COG Functional Category
# see https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2 for other columns that you could use

mkdir /scratch/sogin/tmp.$JOB_ID/eggnog/new_result/
cd /scratch/sogin/tmp.$JOB_ID/eggnog/
emapper_results=$(ls *emapper.annotations)
for i in $emapper_results; do cut -f 1,2,3,7,8,9,10,11,16,21 $i > new_result/${i}_simplified; done

cd /scratch/sogin/tmp.$JOB_ID/eggnog/new_result/
Rscript $path_to_r_parsing_script 

mkdir /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/final/hiqc_bins_test/emapper_test/
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/final/hiqc_bins_test/emapper_test/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```



