library(data.table)
library(tidyverse)

EPIC_EMseq_combined <- readRDS('/home/longtial/EPIC_EMseq_combined_31dec24.rds')

## tms data ----
tms <- EPIC_EMseq_combined_31dec24[,c(1:58)]

tms$avg_meth <- apply(tms[,4:58], 1, function(x) mean(x, na.rm = TRUE))

## epic data ----
epic <- EPIC_EMseq_combined_31dec24[,c(1:3,66:120)]

epic$avg_meth <- apply(epic[,4:58], 1, function(x) mean(x, na.rm = TRUE))

smoothScatter(tms$avg_meth, epic$avg_meth, nrpoints = 1000)
abline(lm(y~x),lty=2)
