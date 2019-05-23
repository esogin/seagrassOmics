# Visualize binning results 
## EM Sogin
## Created May 19,2019
## Updated: May 19, 2019


Description: Visaulize binning results with both GB tools to get cov. by gc stats and differential coverage. Also, experiment with anvio to see wha this tool can help with in terms of analyses. 

## 1. Anvio

```bash
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
```

Notes: 

1. ANvi profiling failed with cpus set to 24, retry this appraoch with cpus set to 4.


## 2. gbtoolscp data

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

## After identification of target bins make a new file called meadow.bins
~/tools/genome-bin-tools-master/accessory_scripts/parse_bin_fasta_files.pl -i meadow.bins -o meadow.bins.tsv 

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

#target metadow bins
d.meadows<-importBins(x=d, file="data/meadow.bins.tsv", to.list=T)


#Plot across libraries

libs<-c('3847_A','3847_B','3847_C','3847_D','3847_E','3847_F','3847_G','3847_H','3847_I') 

bl<-list(d.bins.acido,d.bins.actino, d.bins.alpha, d.bins.alpha.kilo,d.bins.alpha.rhizo,d.bins.alpha.rhodo,d.bins.archaea,d.bins.bact,d.bins.desulfo,d.bins.gamma.gr1,d.bins.gamma.gr2,d.bins.gamma.gr3,d.bins.gamma.gr4,d.bins.gamma.pseudomonadales,d.bins.myo,d.bins.other,d.bins.other2)
lists<-c('acido','actino','alpha','alpha.kilo', 'alpha.rhizo','alpha.rhodo','archaea','bacteroida','desulfobacteria','gamma_gr1','gamma_gr2','gamma_gr3','gamma_gr4','gamma_pseudomonadales', 'myxococcota','other','other2')

# Create plots for all lists --> inspect plots
for (l in 1:length(bl)){
pdf(file=paste("gc_cov_plot_",lists[l],sep='',".pdf"),width=15, height=15)
par(mfrow=c(3,3))
	for(i in 1:9){
		multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=i,cutoff=10000, main=libs[i],legend=T, ylim=c(0.0001,1E3))
	}
	dev.off()
}


#Do differential coverage plots for all bins

for (l in 1:length(bl)){
pdf(file=paste("diff_cov_plot_",lists[l],sep='',".pdf"),width=15, height=15)
par(mfrow=c(3,3))
#Edge slice 4
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(4,7),cutoff=10000,legend=F, ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(4,8),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(4,9),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
#Edge slice 5
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(5,7),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(5,8),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(5,9),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
#Edge Slice 6
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(6,7),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(6,8),cutoff=10000,legend=F,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
multiBinPlot(d, bins=bl[[l]], binNames=names(bl[[l]]),slice=c(6,9),cutoff=10000,legend=T,ylim=c(0.0001, 1E03),xlim=c(0.0001, 1E03))
abline(0,1)
dev.off()
}


# For Target Bins
pdf(file="gc_cov_plot_meadow_bins.pdf",width=12, height=12)
par(mfrow=c(3,3))
multiBinPlot(d, bins=d.meadows, binNames=names(d.meadows),slice=1,cutoff=5000, main=libs[1],legend=T, ylim=c(0.0001,1E3))
for(i in 2:9){
		multiBinPlot(d, bins=d.meadows, binNames=names(d.meadows),slice=i,cutoff=5000, main=libs[i],legend=F, ylim=c(0.0001,1E3))
	}
dev.off()


#
save.image('gbtools.RData')
```

I want to try to make violin plots to show differences in library coverage across habitats

reformat covstats files for inport into R

```bash
libs=$(echo 3847_{A..I}_COVSTATS)
for i in $libs; do cut -f 1,2,3 $i > "$i"_forR;done
```


```R
# Import COV Stats Files 
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
bins<-read.table('data/meadow.bins.tsv', header=F)
colnames(bins)<-c('name','contigID')

#match to names
covstats$name<-bins[match(covstats$ID,bins$contigID),'name'] 
covstats2<-covstats[complete.cases(covstats),]

#Save for plotting
write.csv(covstats2,'results/annotated_covstats_meadow_bins.csv')

```






