#################################
#COMPARING EMseq (TTMS) AND EPIC#
#################################



#########
# work with bsseq object & filter
#########
ml GCC/11.3.0 OpenMPI/4.1.4 R-bundle-Bioconductor/3.15-R-4.2.1

library(comethyl)
library(bsseq)
library(data.table)
data_emseq<-readRDS('/path/to/bsseq/object/from/EMseq/bsseq.rds')
round(colMeans(getCoverage(data_emseq)), 1)
sum(rowSums(getCoverage(data_emseq)) == 0)

options(matrixStats.useNames.NA = "deprecated")
#Filtering CpG-level data for loci with at least X reads in at least XX% of samples
data_filt<-filterCpGs(
  data_emseq,
  cov = "X",
  perSample = 0."X",
  save = FALSE,
  verbose = TRUE)
  #691378 methylation loci
  #55 samples

###
names<-as.data.frame(gsub("/path/to/merged.cov.gz/files/or/CX_report.txt/files/.<file_extension>","",sampleNames(data_filt)))
names(names)[1]<-'names'
setDT(names, keep.rownames = TRUE)
colnames(names) <- c("bsseq_order", "names")
names$bsseq_order = as.numeric(as.character(names$bsseq_order))
names$names <- gsub("XX","X", names$names) 	#to simplify the column names so that they match those in the EPIC data file
names2=read.delim('file_names_file.txt',header=F)
names=merge(names,names2,by.x='names',by.y='V1')
names <- as.data.frame(names)
names=names[order(names$bsseq_order),]

# array data
data=fread('/path/to/EPICdata/output/file/TableControl_SampleMethylationProfile.txt')
probes=read.csv('/path/to/probes/provided/by/EPIC/EPIC-8v2-0_A1_v2.csv')
probes$chr_loc<-paste(probes$CHR,probes$MAPINFO,sep='_')
probes$end<-probes$MAPINFO+1
names(probes)[4]<-'start'
names$array_name<-paste(names$V2,'.AVG_Beta',sep='')

# array sites
array_probes<-probes[which(probes$start>0),c('CHR','start','end')]
array_probes<-array_probes[order(array_probes$CHR,array_probes$start),]
meth<-(getMeth(data_filt, regions = GRanges(array_probes), type = c("raw"),
  what = c("perRegion"), confint = FALSE,
  withDimnames = FALSE))

# bad probes
data_mean<-apply(data[,names$array_name,with=F],1,function(x) mean(x,na.rm=T))
em_mean<-as.data.frame(apply(meth,1,function(x) mean(x,na.rm=T)))
em_mean2<-cbind(em_mean, array_probes)
names(em_mean2)[1]<-'em_mean'
em_mean2$chr_loc<-paste(em_mean2$CHR,em_mean2$start,sep='_')
data_mean2<-as.data.frame(cbind(data_mean,data$TargetID))
names(data_mean2)[2]<-'IlmnID'
data_mean2<-merge(data_mean2,probes[,c('IlmnID','chr_loc')],by='IlmnID')
both<-merge(em_mean2,data_mean2,by='chr_loc')
both<-subset(both,em_mean> -1 & data_mean> -1)
write.table(both,'probe_comparison.txt',row.names=F,sep='\t',quote=F)

# all paired
# variable in both
var<-subset(both,em_mean<0.9 & data_mean<0.9 & em_mean>0.1 & data_mean>0.1)
both=fread('probe_comparison.txt')

# to exclude
exclude <- subset(both, (em_mean>0.75 & data_mean<0.25) | (em_mean<0.25 & data_mean>0.75) )
both<-subset(both, !chr_loc %in% exclude$chr_loc )
# summary(lm(both$em_mean ~ both$data_mean))

# Adjusted R-squared:  XX

# correlation between array and EM-seq
names$beta10x<-NA
names$r.squared10x<-NA
names$beta10x_var<-NA
names$r.squared10x_var<-NA
for (i in 1:dim(names)[1]){
  both1<-merge(data[,c('TargetID',names$array_name[i]),with=F],both[,c('chr_loc','IlmnID')],by.x='TargetID',by.y='IlmnID')
  tmp<-as.data.frame(meth[,i])
  tmp$chr_loc<-paste(array_probes$CHR,array_probes$start,sep='_')
  both1<-merge(both1,tmp,by='chr_loc')
  names$beta10x[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$coefficients[2,1]
  names$r.squared10x[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$r.squared
  both1<-subset(both1,chr_loc %in% var$chr_loc)
  names$beta10x_var[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$coefficients[2,1]
  names$r.squared10x_var[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$r.squared
  print(i)
}
write.table(names,'date_correlation_filter_initials.txt',row.names=F,sep='\t')
pdf(file='date_correlation_filter_initials.pdf')
par(mfrow=c(3,3))
for (i in 1:3){
  both1<-(merge(data[,c('TargetID',names$array_name[i]),with=F],both[,c('chr_loc','IlmnID')],by.x='TargetID',by.y='IlmnID'))
  tmp<-as.data.frame(meth[,i])
  tmp$chr_loc<-paste(array_probes$CHR,array_probes$start,sep='_')
  both1<-as.data.frame(merge(both1,tmp,by='chr_loc'))
  smoothScatter(both1[,3],both1[,4],xlab='EPIC',ylab='EM-seq')
  x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)
  both1<-subset(both1,chr_loc %in% var$chr_loc)
  smoothScatter(both1[,3],both1[,4],xlab='EPIC',ylab='EM-seq')
  x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)
}
dev.off()

# correlation between array and EM-seq - randomize
names$shuf<-sample(names$array_name)
names$beta10x_rand<-NA
names$r.squared10x_rand<-NA
names$beta10x_var_rand<-NA
names$r.squared10x_var_rand<-NA
for (i in 1:dim(names)[1]){
  both1<-merge(data[,c('TargetID',names$shuf[i]),with=F],both[,c('chr_loc','IlmnID')],by.x='TargetID',by.y='IlmnID')
  tmp<-as.data.frame(meth[,i])
  tmp$chr_loc<-paste(array_probes$CHR,array_probes$start,sep='_')
  both1<-merge(both1,tmp,by='chr_loc')
  names$beta10x_rand[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$coefficients[2,1]
  names$r.squared10x_rand[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$r.squared
  both1<-subset(both1,chr_loc %in% var$chr_loc)
  names$beta10x_var_rand[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$coefficients[2,1]
  names$r.squared10x_var_rand[i]<-summary(lm( as.matrix(both1[,3]) ~ as.matrix(both1[,4])))$r.squared
  print(i)
}
write.table(names,'data_correlation_filter_shuffled_initials.txt',row.names=F,sep='\t',quote=F)