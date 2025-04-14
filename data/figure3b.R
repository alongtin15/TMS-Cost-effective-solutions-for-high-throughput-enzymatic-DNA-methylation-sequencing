library(ggplot2)
library(data.table)
library(tidyverse)
library(corrr)

corr <- read_delim('6June24_cai_correlation_5x_ALL.txt')

ggplot(data = corr) +
  geom_density(fill = "#396EC2", alpha = 0.65, aes(x = r.squared10x)) +
  geom_density(fill = "#8ECAFD", alpha = 0.65, aes(x = r.squared10x_var)) +
  theme_classic(base_size = 40) +
  xlab("R-squared") +
  ylab("Density") +
  xlim(0.7,1)
