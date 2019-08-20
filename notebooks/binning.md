# Binning Scripts
## EM Sogin
## May 4 2019
## Updated: Aug 15, 2019

Binning sediment metagenome coassembly using 3 different tools and combining results with DAS tool. 

## 1. Map individual library reads to coassembly
Files are now stored in coassembly folder. Best approach to avoid network congestion: copy all files needed to scratch and then remove prior to rsyncing. Likely a better approach but this works for now to get good stats files. Make sure your coassembly headers are what you want them to look like for all down stream analyses. Suggestion: remove all uncessary info using: 

```bash
sed "s/_flag.*//g" coassembly.contigs.fixed.fa > contigs_use.fa
```
move contigs_use.fa to correct folder for binning and rename as contigs.fa
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# Mapping script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
lib=LIB;
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
cd /scratch/sogin/tmp.$JOB_ID/
#
# Map reads back to coassembly
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz ./
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/assembly/spades/corrected/${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz ./
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/assembly/contigs.fa ./
#
in1=${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz;
in2=${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
#
mkdir /scratch/sogin/tmp.$JOB_ID/bam/
bbmap.sh in=$in1 in2=$in2 ref=contigs.fa covstats="$lib"_COVSTATS out=bam/"$lib".bam scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
samtools sort -o bam/${lib}_sorted.bam bam/"$lib".bam
samtools index bam/${lib}_sorted.bam 
#
#remove extra files
rm contigs.fa ${lib}_highfreq_kmers_1.00.0_0.cor.fastq.gz ${lib}_highfreq_kmers_2.00.0_0.cor.fastq.gz;
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maps/${lib}/;
#rm /scratch/sogin/tmp.$JOB_ID -R;
#
echo "job finished: "
echo "clean up scratch: sogin/tmp.$JOB_ID"
date
```
## 2. Bin coassembly
Using metabat, concoct, maxbin2 for binning and combing bins with DAS tool. 

### Metabat binning
Longest step in process is generating the dpeth file. Move/copy all sorted bam files to metabat folder for binning before running script. 

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
# Binning script with metabat
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/metabat/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID/
# run metabat2
jgi_summarize_bam_contig_depths --minContigLength 1000 --outputDepth depth.txt --pairedContigs paired.txt bam/*sorted.bam 
metabat2 -i contigs.fa -a depth.txt -o result/metabat -v 
cd result/
~/tools/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fa > ../metabat.scaffolds2bins.tsv

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/metabat/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

### Concoct Binning
Copy all BAM files to concoct folder for binning
Need to have installed concoct as a conda environment

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# Binning script with Coconct 1.0.0
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/concoct/ /scratch/sogin/tmp.$JOB_ID; 
cd /scratch/sogin/tmp.$JOB_ID
#
conda activate concoct_env
cut_up_fasta.py contigs.fa  -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed bam/*sorted.bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -t 48 -b result/
merge_cutup_clustering.py result/clustering_gt1000.csv > result/clustering_merged.csv
mkdir result/fasta_bins
extract_fasta_bins.py contigs.fa result/clustering_merged.csv --output_path result/fasta_bins
#
conda deactivate concoct_env
cd result/
sed "s/,/\t/g" clustering_merged.csv > concoct.s2b.tsv
tail -n +2 concoct.s2b.tsv > file.tmp && mv file.tmp concoct.s2b.tsv
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/concoct/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```
###  Maxbin2 binning
This is the longest binning method. It takes 9-10 days to complete. Copy covstats files from the mapping into maxbin2 folder. The covstats have to be reformed for maxbin to just have the contig name and the abundance info. (use ```cut -f 1,5 FILE```)

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# binning script with maxbin2
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maxbin/ /scratch/sogin/tmp.$JOB_ID;
cd /scratch/sogin/tmp.$JOB_ID/
#
# Run maxbin2
run_MaxBin.pl -thread 48  -contig contigs.fa -out maxbin -abund 3847_A_COVSTATS_forMaxBin.txt -abund2 3847_B_COVSTATS_forMaxBin.txt -abund3 3847_C_COVSTATS_forMaxBin.txt -abund4 3847_D_COVSTATS_forMaxBin.txt -abund5 3847_E_COVSTATS_forMaxBin.txt -abund6 3847_F_COVSTATS_forMaxBin.txt -abund7 3847_G_COVSTATS_forMaxBin.txt -abund8 3847_H_COVSTATS_forMaxBin.txt -abund9 3847_I_COVSTATS_forMaxBin.txt
rsync -a /scratch/sogin/tmp.$JOB_ID /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maxbin/
cd result/
~/tools/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fasta > ../maxbin2.scaffolds2bin.tsv
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

## 3. Combing bins using DAS Tool
Prep for das tool
Issue with running das tool: 
The tool works using the test data provided (no issues with installation and dependencies). The tool works iwth a reduced dataset (10K contigs) but fails on whole data. Likely issue: there is an offending sequence in the contigs file. 

My solution: select all headers that were successfully binned across all three approaches. Combing headers into one file and then subset contigs fasta with seqtk to get only the contigs that were successfully binned. Use these to run with DAS_Tool. 
```bash
## select only headers that were binned in assembly 
cut -f 1 concoct.s2b.tsv > fah.concont.txt
cut -f 1 maxbin2.scaffolds2bin.tsv > fah.maxbin.txt
cut -f 1 metabat.s2b.tsv > fah.metabat.txt

cat fah.concont.txt fah.maxbin.txt fah.metabat.txt |sort | uniq > fasta_headers.txt
seqtk subseq contigs.fa fasta_headers.txt > contigs.use.fa 

```

Combine with DAS_Tool
Ran this using a qrsh session with 48 cores, this completed in 1-2 hrs. 
```bash
## Run DAS_Tool
DAS_Tool -i data/concoct.s2b.tsv,data/maxbin2.scaffolds2bin.tsv,data/metabat.s2b.tsv -l concoct,maxbin,metabat -c data/contigs.use.fa  -o result/DASToolRun --threads 48 --search_engine diamond --write_bins 1 
```

## 4. Check quality of bins
Run CheckM, Assign taxonmy to bins using gtdbtk tool kit, run barnap across all bins

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
# run checkm and gtdbtk took kit
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
checkm lineage_wf -t 48 -x fa ./bins checkm_result -f ./checkm.txt
gtdbtk classify_wf --genome_dir ./bins --cpus 48 --out_dir gtdbtk_out -x fa
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

This script runs checkm using more info to get better quality data for the qa pipeline
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
# run checkm
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/checkm/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
checkm tree bins/ checkm_tree/ -t 24 -x .fa
checkm tree_qa ./checkm_tree -o 2 -f tree_qa 
checkm lineage_set ./checkm_tree markers
checkm analyze markers ./bins checkm_analyze -t 24 -x .fa
checkm qa markers ./checkm_analyze -o 2 -t 24 -f checkm_qa_results.txt
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/checkm/ 
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```
Split bins into bac/ and arc/ directories according to taxonomic assignment 
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#
#Run Barrnap 
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/barnap/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
for i in arc/*.fa; do 
	barrnap $i --threads 48 --kingdom arc --outseq ${i%.fa}_rrna.fa;
done;
#
for i in *.fa; do 
	barrnap $i --threads 48 --kingdom bac --outseq ${i%.fa}_rrna.fa;
done;
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/result/barnap/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

Check for tRNAs in finalized bins 

```bash


for n in *.fa; do aragorn -t -fon $n  > trnas/"$n"_trnas; done

grep ">" metabat.6.contigs.fa_trnas | cut -d' ' -f2 |sed 's/([a-z]*)//g'| sed 's/tRNA-//g' | sed 's/?(//g' | sed 's/)//g' | sed 's/|/\n/g' | sort -u | wc -l


# all in one line:
for n in *.fa; do echo $n; aragorn -t -fon $n | grep '>' | cut -d' ' -f2 |sed 's/([a-z]*)//g'| sed 's/tRNA-//g' | sed 's/?(//g' | sed 's/)//g' | sed 's/|/\n/g' | sort -u | wc -l; done > trnas_total


awk 'NR % 2 ==0' trnas_total

```


## END
