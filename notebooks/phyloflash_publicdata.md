# Public datasets

Run phyloflash 

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
