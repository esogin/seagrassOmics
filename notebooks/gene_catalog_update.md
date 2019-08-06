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


##Need to re-visit clustering approach (try number 3)
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
|nr genes/proteins after progidal ||
|nr. proteins after mmseqs (similarity threshold 50% identity)||



Get rep sequence headers and subset protein file 

```bash
grep ">" clusters_rep_seq.fasta > headers
```
* subset protein sequences with list of fasta headers


```bash
seqtk subseq genes/orfs.fa rep_seq_headers > orfs_rep_seqs.fa
```
Nr. seqs: 
17947592

###5 Annotate gene catalog
Using the MPI servers

Try spliting up file 
```bash
split -l 200000 -a 3 -d orfs_rep_seqs.fa input_file.chunk
# generate all the commands that should be distributed in the cluster
for i in *chunk*; do  python2 ./../emapper.py -i $i -m diamond --no_annot --no_file_comments --cpu 24 -o $i; done


cat *.chunk_*.emapper.seed_orthologs > input_file.emapper.seed_orthologs
emapper.py --annotate_hits_table input.emapper.seed_orthologs --no_file_comments -o output_file --cpu 10
```






```bash
python2 ./emapper.py -i my_proteins/orfs_rep_seqs.fa -o catalog_maNOG -m diamond --cpu 24 
```




###4-alt. Annotate gene catalog with emapper.py




-------

*testing*

while running prodigal, run the following tests on subset of gene catalog (old)

4. Construct a NR gene set 
 

* cdhit

```bash
cd-hit-est -i genes_test.fa -o nr_genes.fa -c 0.95 -T 12 -aS 0.9 -n 8 -M 0 -G 0
```
In my test run, took out 5 sequences for the first 16537 genes 
CPU time of 216 (100000)

* remove gene sequences with less then 100 bp

```bash
reformat.sh in=nr_genes.fa minlength=100 out=nr_genes_long.fa
```
removed ~ 6 % of the genes 

* get out fasta headers

```bash
grep ">" nr_genes_long.fa | sed "s/>//g" > headers
```
* subset protein sequences with list of fasta headers

```bash
seqkit grep -n -f headers orfs.faa > orfs_reduced.faa
```

* Annotate with emapper.  
```bash
python2 ./opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper/emapper.py -i orfs_reduced.faa --output orfs_maNOG -m diamond --usemem --cpu 6 ;

python2 /opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper/emapper.py -i orfs_reduced.faa --output orfs_maNOG -m diamond --usemem --cpu 24;


python2.7 /opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper/emapper.py -i test2.fa -m diamond --output test_2 --usemem --cpu 24 

python2.7 ./emapper.py -i my_proteins/test.fa -o test4 -m diamond --data_dir ./data/

python2 /opt/extern/bremen/symbiosis/sogin/tools/eggnog-mapper/emapper.py -i metabat.98.contigs_orfs.faa --output test_maNOG -m diamond --usemem --cpu 24;
```




