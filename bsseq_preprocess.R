##########################
#PREPROCESSING BSSEQ IN R#
##########################

library(tidyverse)
library(data.table)
library(bsseq)

bsseq_RDS <- readRDS('/path/to/generated/bsseq/object/bsseq.rds')

## getting the data frames workable ----
#grabbing the chromosomes and the positions
chrom <- seqnames(rowRanges(bsseq_RDS))
pos <- start(rowRanges(bsseq_RDS))

#getting the methylation and coverage files
methylation <- as.data.frame(getMeth(bsseq_RDS, type = "raw"))
coverage <- as.data.frame(getCoverage(bsseq_RDS, type = "Cov"))

#merging the dataframes
methylation <- data.frame(chromosome = chrom, position = pos, methylation)
coverage <- data.frame(chromosome = chrom, position = pos, coverage)

# filtering for coverage
library(comethyl)
filterCpG(bsseq_RDS,
	cov = "X",
	perSample = "X") #then would have to process with the code from above to get the coverage and methylation matrices

	#OR#

filt_df <- coverage[apply(coverage[,"X":"XX"], 1, function(x){(sum(x>"X", na.rm = T))/ncol(coverage)}) >=0."XX",]
	#replace the "X"'s with the specific rows with the coverage information, the coverage value you are filtering for, and the percentage of samples, respectively