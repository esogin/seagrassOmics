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





