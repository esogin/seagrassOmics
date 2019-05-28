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
			phyloFlash.pl -lib ${i%%.fna} -CPUs 1 -read1 $i;
		rm core
	done		 
echo "job finished: "
date
```

```bash
values=$(ls *fna);
for i in $values;do nohup phyloFlash.pl -lib ${i%%.3.fna} -CPUs 1 -read1 fna/$i -zip > phyloflash.log; done	
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




