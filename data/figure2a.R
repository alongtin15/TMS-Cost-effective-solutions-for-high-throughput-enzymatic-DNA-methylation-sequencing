library(rlang)
library(ggplot2)
library(tidyverse)

metadata <- read.csv('~/Documents/vandy/EMseq/DNA methylation metadata - plexing.csv')
metadata$Mapping_efficency <- as.numeric(gsub("[\\%,]", "", metadata$Mapping_efficency))
transform(metadata, Mapping_efficency = as.numeric(Mapping_efficency))
metadata$methylation_CHH <- as.numeric(gsub("[\\%,]", "", metadata$methylation_CHH))
transform(metadata, methylation_CHH = as.numeric(methylation_CHH))

ggplot(data = metadata, aes(x=as.factor(Plexing), y=Mapping_efficency, group = as.factor(Plexing))) +
  geom_violin(alpha = 0.45, fill = "#123067") +
  geom_boxplot(alpha = 0.65, fill = "#123067", width = 0.15) +
  geom_jitter(alpha = 0.75, color = "#123067", size = 6) +
  ylab('Mapping Efficiency (%)') +
  xlab('Plexing Type') +
  #labs(title = 'Mapping Efficiency by Plexing Type') +
  theme(legend.position = "none") +
  theme_classic(base_size = 40) +
  ylim(40,90)
