# Public datasets

Run phyloflash on other data (e.g., mgrast)
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -V
#$ -q main.q@@himem
#
# Assembly script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
cd /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/soil_metagenomes/mgrast_datasets/fraser/
values=$(ls *fna)
for i in $values; 
	do 
		echo $i 
			phyloFlash.pl  - lib ${i%%.fa} -CPUs 24 -read1 $i;
		rm core
	done		 
echo "job finished: "
date
```

Run phyloflash on ENA Data

```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 64
#$ -V
#$ -q main.q@@himem
#
# Assembly script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
cd /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/soil_metagenomes/public_data/
value=$(<accession_numbers)

for i in $values; 
	do 
		echo $i 
			ENA_phyloFlash.pl --acc $i  --cleanup --phyloFlash "taxlevel 3";
			rm core
	done		 
echo "job finished: "
date
```
Compare phyloflash results 
```bash
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 64
#$ -V
#$ -q main.q@@himem
#
# Assembly script
echo "job started: " 
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
cd /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/soil_metagenomes/public_data/done
phyloFlash_compare.pl --allzip --out public_data --keeptmp --task barplot
echo "job finished: "
date
```




