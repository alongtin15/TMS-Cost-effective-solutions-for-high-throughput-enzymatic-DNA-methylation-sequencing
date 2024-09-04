#############################
#R COMMANDS FOR MAKING BSSEQ#
#############################

#place this in a file in your directory (ie. bsseq.R)

library(Biostrings)
library(bsseq)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

chrs <- names(Hsapiens)[1:22]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))

twist_regions<-read.delim("/path/to/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed",header = F)
twist_regions<-twist_regions[,1:3]
colnames(twist_regions) <- c("chr","start","end")
twist_regions2<-twist_regions
twist_regions2$start<-twist_regions$start-200
twist_regions2$end<-twist_regions$end+200

sites_in_twist <- as.data.frame(GenomicRanges::intersect(x = cpgr,
                                                         y = GRanges(twist_regions2)))
sites_in_twist2<-subset(sites_in_twist,width==2)
sites_in_twist2$end<-sites_in_twist2$start

sites_in_twist3<-sites_in_twist2
sites_in_twist3$start<-sites_in_twist3$end<-sites_in_twist2$end+1

sites_in_twist2$strand<-'+'
sites_in_twist3$strand<-'-'
sites_in_twist<-rbind(sites_in_twist2,sites_in_twist3)
sites_in_twist<-sites_in_twist[order(sites_in_twist$seqnames,sites_in_twist$start),]
sites_in_twist$width<-1

bismarkBSseq <- read.bismark(files = list.files("/path/to/CG/files/",
                                                full.names = T, pattern=".CG_report.txt"), #".CG_report.txt" can be replaced with .cov.gz#
                             strandCollapse = TRUE, loci=GRanges(sites_in_twist),
                             verbose = TRUE, BACKEND = "HDF5Array",
                             dir="/path/to/output/directory/bsseq/", replace=TRUE, nThread=4)

saveRDS(bismarkBSseq,file="bsseq_object.rds")
q(save="no")