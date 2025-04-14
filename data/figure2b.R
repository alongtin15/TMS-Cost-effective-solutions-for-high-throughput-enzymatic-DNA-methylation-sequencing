library(rlang)
library(ggplot2)
library(tidyverse)

metadata <- read.csv('DNA methylation metadata - input.csv')
metadata$Mapping_efficency <- as.numeric(gsub("[\\%,]", "", metadata$Mapping_efficency))
transform(metadata, Mapping_efficency = as.numeric(Mapping_efficency))
metadata$methylation_CHH <- as.numeric(gsub("[\\%,]", "", metadata$methylation_CHH))
transform(metadata, methylation_CHH = as.numeric(methylation_CHH))

##mapping efficiency x input type (boxplot)
ggplot(data = metadata, aes(x=reorder(Input, Mapping_efficency), y=Mapping_efficency, group = Input)) +
  geom_jitter(alpha = 0.85, color = "#447121", size = 10, position = position_jitter(width = 0.0, height = 0.1)) +
  ylab('Mapping Efficiency (%)') +
  xlab('Input Amount (ng)') +
  #labs(title = 'Mapping Efficiency by Input Amount') +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  scale_fill_discrete(limits = c("25","50","100","200","400")) +
  theme_classic(base_size = 40) +
  ylim(40,90)
