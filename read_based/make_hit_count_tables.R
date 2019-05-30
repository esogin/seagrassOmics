#!/usr/bin/Rscript

#Cherries
files<-list.files(pattern="matches_cherry.m8", recursive = T)
files<-files[grep('function_diamond',files)]
fnames<-gsub(pattern=".*function_diamond/","",x=files)
df<-read.table(files[1],sep="\t", header=F)
colnames(df)<-c('qseqid', 'sseqid', 'pident' ,'length', 'mismatch', 'gapopen' ,'qstart', 'qend' ,'sstart' ,'send' ,'evalue', 'bitscore')
counts<-data.frame(table(df$sseqid))
colnames(counts)<-c('Hit',paste0(fnames[1],'_frq',sep=""))
for (i in 2:length(files)){
df2<-read.table(files[i],sep="\t", header=F)
colnames(df2)<-c('qseqid', 'sseqid', 'pident' ,'length', 'mismatch', 'gapopen' ,'qstart', 'qend' ,'sstart' ,'send' ,'evalue', 'bitscore')
counts2<-data.frame(table(df2$sseqid))
colnames(counts2)<-c('Hit',paste0(fnames[i],'_frq',sep=""))
counts<-merge(counts, counts2, by="Hit", all=T)
}
counts[is.na(counts)]<-0
write.csv(x=counts,'results/read_based/hit_counts_cherries.csv')

#cazy
files<-list.files(pattern="matches_cazy.m8", recursive = T)
files<-files[grep('function_diamond',files)]
fnames<-gsub(pattern=".*function_diamond/","",x=files)
df<-read.table(files[1],sep="\t", header=F)
colnames(df)<-c('qseqid', 'sseqid', 'pident' ,'length', 'mismatch', 'gapopen' ,'qstart', 'qend' ,'sstart' ,'send' ,'evalue', 'bitscore')
counts<-data.frame(table(df$sseqid))
colnames(counts)<-c('Hit',paste0(fnames[1],'_frq',sep=""))
for (i in 2:length(files)){
df2<-read.table(files[i],sep="\t", header=F)
colnames(df2)<-c('qseqid', 'sseqid', 'pident' ,'length', 'mismatch', 'gapopen' ,'qstart', 'qend' ,'sstart' ,'send' ,'evalue', 'bitscore')
counts2<-data.frame(table(df2$sseqid))
colnames(counts2)<-c('Hit',paste0(fnames[i],'_frq',sep=""))
counts<-merge(counts, counts2, by="Hit", all=T)
}
counts[is.na(counts)]<-0
write.csv(x=counts,'results/read_based/hit_counts_cazy.csv')

#sugars
files<-list.files(pattern="matches_sugars.m8", recursive = T)
files<-files[grep('function_diamond',files)]
fnames<-gsub(pattern=".*function_diamond/","",x=files)
df<-read.table(files[1],sep="\t", header=F)
colnames(df)<-c('qseqid', 'sseqid', 'pident' ,'length', 'mismatch', 'gapopen' ,'qstart', 'qend' ,'sstart' ,'send' ,'evalue', 'bitscore')
counts<-data.frame(table(df$sseqid))
colnames(counts)<-c('Hit',paste0(fnames[1],'_frq',sep=""))
for (i in 2:length(files)){
df2<-read.table(files[i],sep="\t", header=F)
colnames(df2)<-c('qseqid', 'sseqid', 'pident' ,'length', 'mismatch', 'gapopen' ,'qstart', 'qend' ,'sstart' ,'send' ,'evalue', 'bitscore')
counts2<-data.frame(table(df2$sseqid))
colnames(counts2)<-c('Hit',paste0(fnames[i],'_frq',sep=""))
counts<-merge(counts, counts2, by="Hit", all=T)
}
counts[is.na(counts)]<-0
write.csv(x=counts,'results/read_based/hit_counts_sugars.csv')
