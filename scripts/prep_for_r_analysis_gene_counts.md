
#Prepare scafstats files for import into R
for i in $scafstats; do sed "s/ # /_/g" $i > fixed_"$i"; done



```r
## Combine read counts into single file for profiles
files<-list.files(pattern="fixed.*.txt")

df<-read.table(files[1],sep="\t", header=F)
colnames(df)<-c('name','%unambiguousReads',	'unambiguousMB','%ambiguousReads','ambiguousMB','unambiguousReads','ambiguousReads',	'assignedReads','assignedBases')
count<-data.frame(df$name, as.numeric(df$assignedReads))

colnames(count)<-c('name',paste0(files[1],sep=""))
#count$assignedReads<-as.vector(count$assignedReads)

for (i in 2:length(files)){
df2<-read.table(files[i],sep="\t", header=F)
colnames(df2)<-c('name','%unambiguousReads','unambiguousMB','%ambiguousReads','ambiguousMB','unambiguousReads','ambiguousReads',	'assignedReads','assignedBases')
count2<-data.frame(df2$name, as.numeric(df2$assignedReads))
colnames(count2)<-c('name',paste0(files[i],'_count',sep=""))
count<-merge(count, count2, by="name", all=T)
}

count<-count[is.na(count)]<-0
count$Total<-rowSums(count[,grep("fixed",colnames(count))])
range(count$Total)
summary(count$Total)

colSums(count[,grep("fixed",colnames(count))]) #between 49,820,708 and 101,056,230 reads maped to catalog

count.abund<-count[count$Total > 10, ]


## Add in gene identifiers to profile table
annotations<-read.table('~/extern/Data/SedimentMG/processed_reads/libraries/library_3847/gene_catalog2/annotate_me/output_orfs/h.out_modified.csv',sep=',',header=T)
colnames(annotations)[1]<-'domain'

head(count.abund)
colnames(count.abund)[1]<-'name'
count.abund$name<-as.vector(count.abund$name)
count.abund$name2<-paste(sapply(strsplit(count.abund$name,"_"),'[',1),sapply(strsplit(count.abund$name,"_"),'[',2),sep="_")
count.abund$DomainID<-annotations[match(count.abund$name2,annotations$query_name),'domain']

# write to output
write.csv(count.abund,'gene_counts_cazyAnnot_May52019.csv')
```