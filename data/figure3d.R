## loading libraries and data ----
library(tidyverse)
library(data.table)
library(dplyr)
library(bsseq)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(GenomicRanges)
library(BiocManager)

wgbs_avgmeth <- read.delim('/nobackup/lea_lab/longtial/EMseq23-24/wgbs/wgbs_avgmeth.txt')
EPIC_EMseq_combined_31dec24 <- readRDS("/Users/longtial/Downloads/EPIC_EMseq_combined_31dec24.rds")

##tms data ----
tms <- EPIC_EMseq_combined_31dec24[,c(1:58)]

tms$avg_meth <- apply(tms[,4:58], 1, function(x) mean(x, na.rm = TRUE))

## epic data ----
epic <- EPIC_EMseq_combined_31dec24[,c(1:3,66:120)]

epic$avg_meth <- apply(epic[,4:58], 1, function(x) mean(x, na.rm = TRUE))

## plotting ----
ggplot(data = NULL) +
  geom_density(fill = "#DE3C37", alpha = 0.5, data = wgbs_avgmeth, aes(x = avg_meth)) +
  geom_density(fill = "#123067", alpha = 0.5, data = tms, aes(x = avg_meth)) +
  geom_density(fill = "#8ECAFD", alpha = 0.5, data = epic, aes(x = avg_meth)) +
  theme_classic(base_size = 40) +
  xlab("Average Methylation (%)") +
  ylab("Density")
