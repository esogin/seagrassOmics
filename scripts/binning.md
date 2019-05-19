#Binning Workflow
#Esogin
#May 4 2019
#Updated: May 15, 2019

First binning approach there was an issue with the mapping files. This needs to be fixed to ensure approapriate binning. These scripts now use better mapped data for coverage stats

1. Map library reads to coassembly & generate a depth file for metabat and concoct binning approaches
Done May 15 2019
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
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/congits.fa ./
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
2. Binning approaches for coassembly

# Do metabat binning 
Done May 15 2019
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

# Do concoct
Done May 15 2019
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
# Do maxbin2
Done May 15 2019
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

Prep for das tool
```bash
## select only headers that were binned in assembly 
cut -f 1 concoct.s2b.tsv > fah.concont.txt
cut -f 1 maxbin2.scaffolds2bin.tsv > fah.maxbin.txt
cut -f 1 metabat.s2b.tsv > fah.metabat.txt

cat fah.concont.txt fah.maxbin.txt fah.metabat.txt |sort | uniq > fasta_headers.txt
seqtk subseq contigs.fa fasta_headers.txt > contigs.use.fa 

```


3. Combine with DAS_Tool
Done May 15 2019
```bash
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
# combine bins with dastool
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/ /scratch/sogin/tmp.$JOB_ID;
cd /scratch/sogin/tmp.$JOB_ID/
#
## Run DAS_Tool
DAS_Tool -i data/concoct.s2b.tsv,data/maxbin2.scaffolds2bin.tsv,data/metabat.s2b.tsv -l concoct,maxbin,metabat -c data/contigs.use.fa  -o result/DASToolRun --threads 48 --search_engine diamond --write_bins 1 
#
rsync -a /scratch/sogin/tmp.$JOB_ID /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/dastool/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```
Run CheckM, Assign taxonmy to bins using gtdbtk tool kit, run barnap across all bins
Done May 15 2019
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/drep/drep_bins/dereplicated_genomes/ /scratch/sogin/tmp.$JOB_ID/; 
cd /scratch/sogin/tmp.$JOB_ID/
#
checkm lineage_wf -t 48 -x fa ./bins checkm_result -f ./checkm.txt
checkm qa
gtdbtk classify_wf --genome_dir ./bins --cpus 48 --out_dir gtdbtk_out -x fa
#
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/drep/drep_bins/dereplicated_genomes/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```
Split bins into bac/ and arc/ directories according to taxonomic assignment 
Done May 15 2019
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
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/combined_result/bins/ /scratch/sogin/tmp.$JOB_ID/; 
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
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/bins/combined_result/bins/
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```



4. Visualize binning results 

### 4a With anvio
```bash
conda activate anvio5

#Reformat fasta file of contigs and remove contigs with length less then 1000
#Keep bin names the same so that we can import other mapping data onto contigs database at end
anvi-script-reformat-fasta contigs.fa -o contigs-fixed.fa -l 1000 
mv contigs-fixed.fa contigs.fa
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n "Sediment MG Database"
anvi-run-hmms -c contigs.db --num-threads 48
#anvi-display-contigs-stats contigs.db  # only run if have ssh tunnel in place 
#anvi-setup-ncbi-cogs #do only once
#Annotate assembly with cog categories, may be useful later on
anvi-run-ncbi-cogs -c contigs.db --num-threads 48

#Initialize BAM files
bams=$(echo 3847_{A..I}_sorted.bam)
for sample in $bams; do anvi-init-bam bam/$sample -o $sample; done

#Profile BAM files
bams=$(echo 3847_{A..I})
for i in $bams;do
	anvi-profile -i data/${i}_sorted.bam -c contigs.db --output-dir ${i}_profile --sample-name profile_"$i";
done
```



### 4b with gbtools
prep covstats files for gbtools
prep of bin files
```bash

~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i archaea.bins -o archaea.bins.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i acidobacteria.bins -o acidobacteria.bins.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i actinobacteria.bins -o actinobacteria.bins.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i bacteroida.bins -o bacteroida.bins.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i desulfobacteria.bins -o desulfobacteria.bins.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i myxococcota.bins -o myxococcota.bins.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i alpha.kilo -o alpha.kilo.table.tsv
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i alpha.rhizobales.bins -o alpha.rhizo.bins.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i alpha.rhodo -o alpha.rhodo.bins.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i alpha.bins -o alpha.bins.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i gamma.gr1.bins -o gamma.gr1.bins.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i gamma.gr2.bins -o gamma.gr2.bins.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i Gamma.Pseudomonadales -o Gamma.Pseudomonadales.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i gamma.gr3 -o gamma.gr3.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i gamma.gr4 -o gamma.gr4.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i other.bins -o other.bins.tsv 
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i other2.bins -o other2.bins.tsv 

```
```R
library(gbtools)
files<-list.files(pattern="*COVSTATS", recursive=T, include.dirs=T)
#Import covstats files 
d <- gbt(covstats=c("covstats/3847_A_COVSTATS", "covstats/3847_B_COVSTATS", "covstats/3847_C_COVSTATS","covstats/3847_D_COVSTATS","covstats/3847_E_COVSTATS", "covstats/3847_F_COVSTATS","covstats/3847_G_COVSTATS", "covstats/3847_H_COVSTATS","covstats/3847_I_COVSTATS"))
summary(d)

## Import Bins (seperated by phylum and sometimes subphylum level)
d.bins.acido<-importBins(x=d, file="data/acidobacteria.bins.table.tsv", to.list=T)
d.bins.actino<-importBins(x=d, file="data/actinobacteria.bins.table.tsv", to.list=T)
d.bins.alpha<-importBins(x=d, file="data/alpha.bins.tsv", to.list=T)
d.bins.alpha.kilo<-importBins(x=d, file="data/alpha.kilo.table.tsv", to.list=T)
d.bins.alpha.rhizo<-importBins(x=d, file="data/alpha.rhizo.bins.tsv", to.list=T)
d.bins.alpha.rhodo<-importBins(x=d, file="data/alpha.rhodo.bins.tsv", to.list=T)
d.bins.archaea<-importBins(x=d, file="data/archaea.bins.table.tsv", to.list=T)
d.bins.bact<-importBins(x=d, file="data/bacteroida.bins.table.tsv", to.list=T)
d.bins.desulfo<-importBins(x=d, file="data/desulfobacteria.bins.table.tsv", to.list=T)
d.bins.gamma.gr1<-importBins(x=d, file="data/gamma.gr1.bins.tsv", to.list=T)
d.bins.gamma.gr2<-importBins(x=d, file="data/gamma.gr2.bins.tsv", to.list=T)
d.bins.gamma.gr3<-importBins(x=d, file="data/gamma.gr3.tsv", to.list=T)
d.bins.gamma.gr4<-importBins(x=d, file="data/gamma.gr4.tsv", to.list=T)
d.bins.gamma.pseudomonadales<-importBins(x=d, file="data/Gamma.Pseudomonadales.tsv", to.list=T)
d.bins.myo<-importBins(x=d, file="data/myxococcota.bins.table.tsv", to.list=T)
d.bins.other<-importBins(x=d, file="data/other.bins.tsv", to.list=T)
d.bins.other2<-importBins(x=d, file="data/other2.bins.tsv", to.list=T)

#Plot across libraries

libs<-c('3847_A','3847_B','3847_C','3847_D','3847_E','3847_F','3847_G','3847_H','3847_I') 

bl<-list(d.bins.acido,d.bins.actino, d.bins.alpha, d.bins.alpha.kilo,d.bins.alpha.rhizo,d.bins.alpha.rhodo,d.bins.archaea,d.bins.bact,d.bins.desulfo,d.bins.gamma.gr1,d.bins.gamma.gr2,d.bins.gamma.gr3,d.bins.gamma.gr4,d.bins.gamma.pseudomonadales,d.bins.myo,d.bins.other,d.bins.other2)
lists<-c('acido','actino','alpha','alpha.kilo', 'alpha.rhizo','alpha.rhodo','archaea','bacteroida','desulfobacteria','gamma_gr1','gamma_gr2','gamma_gr3','gamma_gr4','gamma_pseudomonadales', 'myxococcota','other','other2')

# Create plots for all lists --> inspect plots
for (l in 1:length(bl)){
pdf(file=paste("gc_cov_plot_",lists[l],sep='',".pdf"),width=15, height=15)
par(mfrow=c(3,3))
	for(i in 1:9){
		multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=i,cutoff=10000, main=libs[i],legend=T)
	}
	dev.off()
}


#Do differential coverage plots for all bins

for (l in 1:length(bl)){
pdf(file=paste("diff_cov_plot_",lists[l],sep='',".pdf"),width=15, height=15)
par(mfrow=c(3,3))
#Edge slice 4
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(4,7),cutoff=10000,legend=F)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(4,8),cutoff=10000,legend=F)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(4,9),cutoff=10000,legend=F)
#Edge slice 5
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(5,7),cutoff=10000,legend=F)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(5,8),cutoff=10000,legend=F)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(5,9),cutoff=10000,legend=F)
#Edge Slice 6
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(6,7),cutoff=10000,legend=F)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(6,8),cutoff=10000,legend=F)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(6,9),cutoff=10000,legend=F)
dev.off()
}

#
save.image('gbtools.RData')
```

