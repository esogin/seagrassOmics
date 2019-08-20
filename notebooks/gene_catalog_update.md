#Gene Catalog Update
## EM Sogin 
## Created July 23, 2019
## Updated August 5, 2019


Description: in our initial gene catalog appraoch, we have issues with annotation and mapping results. This is likely because our inital gene catalog has over 83M genes to filter through. Try to reduce number of reads by filtering out short sequences and reimplement a gene catalog annotation strategy that is more tractable. 

initial commands to dereplicate and cluster sequences
```bash
vsearch --derep_fulllength allcontigs_renamed.fa --output clusters/rep_set.fa  --relabel group --relabel_keep
mmseqs easy-linclust clusters/rep_set.fa clusters tmp
```


###1. Investgiate old gene catalog 

Run quast to get an idea of average sequence lenght and other stats from initial 83M gene catalog. 

```bash
quast.py -f clusters_rep_seq.fasta -o quast_run_repset --fast 
```

###2. Remove short sequences from gene catalog
length of bacteria genes > 500 bp 
http://bioscience.jbpub.com/cells/MBIO137.aspx
as in https://www.nature.com/articles/nmicrobiol2016161#methods

Also, fix reads by adding unscore to name
```bash
reformat.sh in=clusters_rep_seq.fasta minlength=500 out=clusters_rep_500bp.fasta
reformat.sh in=clusters_rep_seq_500bp.fasta out=clusters_rep_seq_500bp_fixed.fasta addunderscore
sed 's/__//g' clusters_rep_seq_500bp_fixed.fasta > clusters_rep_seq_500bp_fixed2.fasta
```

number of reads that went in: 73,357,971
number of reads that came out: 12,599,825

###3. run prodigal 

```blast
prodigal -i clusters_rep_seq_500bp.fasta -o genes/coords.gff  -a genes/orfs.faa  -p meta -d genes/genes.fa -f gff; 
```

number of proteins: 
21,077,130 (25% of original protein file!)

###4. Re-cluster to construct non-redundant gene set

I am not convinced this will reduce the gene catalog in a meaningful way as we already did this before calling the genes. However, proceed to impliement cd-hit and see if we can get a reduction. 

**Cavet: all published gene catalogues do this step here**
either with BLAT or cd-hit. DHI will take a while...

```bash
cd-hit-est -i genes/genes.fa -o genes/nr_genes.fa -c 0.95 -T 12 -aS 0.9 -n 8 -M 0 -G 0;
```


--> cd-hit slow, try working with blat again 

```bash
blat genes.fa genes.fa -minIdentity=95 -t=dna -tileSize=11 -stepSize=10 output.psl
```

--> Blat ran out of memory on himem machine. 


Try reclusting a second time with mmseqs to get non-redundant genes out.

```bash
mmseqs easy-linclust genes/genes.fa clusters tmp
```
--> mmseqs worked very quickly (< 2 mins) 

Nr. sequences in: 
21,077,130

Nr of sequences out: 
17,947,592


##Need to re-visit clustering approach (3rd time is the charm)
This is relatively small reduction in sequences. 

**Update: The workflow for mmseqs was incorrect applied above, re-implement here:**

This workflow takes combined assembilies dereplicates them, reduces the size and then calls proteins and genes on entire catalog. The piepline will then use mmseqs to create a non-redundant gene catalog for analysis. 

my qsub script for a better pipeline

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 48
#$ -V
#$ -q main.q
#$ -l "h=gc-node-22"
#Creating a gene catalog
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
mkdir /scratch/sogin/tmp.$JOB_ID -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/ /scratch/sogin/tmp.$JOB_ID/

cd /scratch/sogin/tmp.$JOB_ID/


echo "reduce catalog to only include sequences with at least 500 bp"
echo "already done"
#reformat.sh in=allcontigs_renamed.fa minlength=500 out=contigs_500bp.fasta

echo "fix fasta headers"
echo "already done"
#reformat.sh in=contigs_500bp.fasta out=contigs_500bp_f1.fasta addunderscore
#sed 's/_flag.*/g' contigs_500bp_f1.fasta > contigs.fasta

echo "dereplicate catalog and run initial clustering using mmseqs"
mkdir /scratch/sogin/tmp.$JOB_ID/clusters/

#vsearch --derep_fulllength contigs.fasta --output clusters/rep_set.fa  --relabel group --relabel_keep

#this clusters to 100 % identity on gene level (redundent step to above - see how many more sequences it removes)
#mmseqs easy-linclust clusters/rep_set.fa clusters/derep_seqs tmp

echo "call genes"
mkdir /scratch/sogin/tmp.$JOB_ID/prodigal/
mkdir /scratch/sogin/tmp.$JOB_ID/prodigal/genes/

prodigal -i clusters/derep_seqs_rep_seq.fasta -o prodigal/genes/coords.gff  -a prodigal/genes/orfs.faa  -p meta -d prodigal/genes/genes.fa -f gff;

mkdir prodigal/mmseqs/
cd prodigal/mmseqs/

mmseqs easy-cluster ../genes/orfs.faa orfs_clusters_50 tmp --min-seq-id 0.5
mmseqs easy-cluster ../genes/orfs.faa orfs_clusters_90 tmp --min-seq-id 0.9

echo "rsync back"
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/

echo "job finished: "
date
echo "CLEAN UP SCRATCH!!!"

if [ "$?" -eq "0" ]; then
  rm /scratch/tmp.$JOB_ID -R;
  echo "Done"
else
  echo "Error while running rsync"
fi
```

Clustering results: 
| Step | number of sequences|
|:-------:|:-------------:|
| nr. contigs (initial)|75,805,676|
| nr. contigs after removing short seqs|13,146,164|
| nr. contigs after derep (vsearch) |13,132,325|
| nr. contigs after mmseqs|12,701,223|
|nr genes/proteins after progidal |21,199,115|
|nr. proteins after mmseqs (similarity threshold 50% identity)|8,930,797|
|nr. proteins after mmseqs (similarity threshold 90% identity)|17,288,787|

*Annotation pipeline of 50% ID catalog*
for some reason I can't get emapper to work on the server. Instead run on MPI servers 

```bash
#
split -l 200000 -a 3 -d orfs_clusters_50_rep_seq.fasta input_file.chunk
# generate all the commands that should be distributed in the cluster

e_path=~/devtools/eggnog-mapper/

for i in *chunk*; do  python2 $e_path/emapper.py -i $i -m diamond --no_annot --no_file_comments --cpu 48 -o $i; done


cat *.chunk*.emapper.seed_orthologs > input_file.emapper.seed_orthologs
python2 $e_path/emapper.py --annotate_hits_table input_file.emapper.seed_orthologs --no_file_comments -o output_file --cpu 48
```

Get rep sequence headers and subset protein file 
using cluster 50 file
```bash
grep ">" orfs_clusters_50_rep_seq.fasta > headers
```
subset gene sequences with list of fasta headers from protein file

```bash
sed 's/>//g' headers > headers_list
seqtk subseq genes.fa headers_list > genes_50_rep_seqs.fa
```
Check if have same number of sequences

Nr. seqs: 
8930797


## 5. Map libraries back to catalog 
This is not as trivial as I had initial hoped. 

Issues: 

1) only 25% of reads map back to catalog

- try to map read files (1 file) to both the 50% clustered catalog and the 90% clustered catalog. If diffent read mapping, this indiciates a sensitivity issue in the mapping results. 


This script will try multiple mapping appraoches including the one described above. Also interested to see if the mapping sucess varies depending on mappers. 
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
#Template script to map libraries to gene catalog
#notes: the libraries being mapped have not been error corrected with spades and are NOT kmer normalized, hence the nn director name (for no normalization)
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/ref_genes/prodigal/mmseqs/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
echo "map to catalog";
#
lib=3847_A;
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/qc/"$lib"_q2_ktrimmed.fq.gz ./;
reformat.sh in="$lib"_q2_ktrimmed.fq.gz out=test_reads.fq.gz reads=1000000

mkdir /scratch/sogin/tmp.$JOB_ID/bbmap/
cd /scratch/sogin/tmp.$JOB_ID/bbmap/

bbmap.sh in=../test_reads.fq.gz ref=../genes_50_rep_seqs.fa covstats=test_covstats_50 scafstats=test_scafstats_50.txt statsfile=test_stderr_50;
bbmap.sh in=../test_reads.fq.gz ref=../genes_90_rep_seqs.fa covstats=test_covstats_90 scafstats=test_scafstats_90.txt statsfile=test_stderr_90;

mkdir /scratch/sogin/tmp.$JOB_ID/bowtie/
cd /scratch/sogin/tmp.$JOB_ID/bowtie/

bowtie2-build ../genes_50_rep_seqs.fa genes_50_index
bowtie2 -x genes_50_index -U ../test_reads.fq.gz -S test_50.sam --sensitive 
samtools view -bS -o test_50.bam test_50.sam
samtools sort -o test_50_sorted.bam test_50.bam
samtools view -F 0x4 test_50_sorted.bam | cut -f 1 | sort | uniq | wc -l 

mkdir /scratch/sogin/tmp.$JOB_ID/bwa/
cd /scratch/sogin/tmp.$JOB_ID/bwa/

bwa index -a bwtsw ../genes_50_rep_seqs.fa 
bwa mem ../genes_50_rep_seqs.fa ../test_reads.fq.gz > test_50_bwa.sam
samtools view -bS -o test_50_bwa.bam test_50_bwa.sam
samtools sort -o test_50_sorted_bwa.bam test_50_bwa.bam
samtools view -F 0x4 test_50_sorted_bwa.bam | cut -f 1 | sort | uniq | wc -l 

#
mkdir /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/mapping_test/
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/"$lib"/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```



More sensitive approaches inclued: 

- Translated read mapping with DiamondX
- hmmer searches of translated read profiles against gene catalog 

If using hmmer, we have the hmmprofiles as a result of mmseqs

```bash
hmmsearch 
```





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
#Template script to map libraries to gene catalog
#notes: the libraries being mapped have not been error corrected with spades and are NOT kmer normalized, hence the nn director name (for no normalization)
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
#
mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/ref_genes/ /scratch/sogin/tmp.$JOB_ID/
cd /scratch/sogin/tmp.$JOB_ID/
#
echo "map to catalog";
#
lib=LIB;
cp /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/$lib/qc/"$lib"_q2_ktrimmed.fq.gz ./;
bbmap.sh in=${lib}_q2_ktrimmed.fq.gz ref=genes_50_rep_seqs.fa covstats="$lib"_covstats out="$lib".bam scafstats="$lib"_scafstats.txt statsfile="$lib"_stderr;
rm *fq.gz
rm genes.fa
rm -r ref/
#
mkdir /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/"$lib"/
rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/"$lib"/;
rm /scratch/sogin/tmp.$JOB_ID -R;
echo "job finished: "
date
```

apply to all libraries to get individual scripts

```bash
libs=$(echo 3847_{A..I})
for i in $libs; 
do 
	cat template_map_to_catalog.sh | sed "s/LIB/$i/g" > map_to_catalog_"$i".sh; 
done

```
Run Feature counts

need to generate an SAF file


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

mkdir /scratch/sogin/tmp.$JOB_ID/ -p; 
rsync -a /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/feature_counts/ /scratch/sogin/tmp.$JOB_ID/

cd /scratch/sogin/tmp.$JOB_ID/

featureCounts -F 'GTF' -t 'CDS' -g 'ID' -T 12 -a coords_fixed.gff -o result 3847_A.bam

rsync -a /scratch/sogin/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog/catalog/maps/feature_counts;
rm /scratch/sogin/tmp.$JOB_ID -R;

echo "job finished: "
date
```



