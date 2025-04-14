library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggforce)

capuchin <- read_csv('NHP metadata - capuchin.csv')
macaque <- read_csv('NHP metadata - macaque.csv')
gelada <- read_csv('NHP metadata - gelada.csv')

species <- "capuchin"
capuchin$species <- species

combined <- rbind(capuchin, macaque, gelada)

combined$mapping_efficiency <- as.numeric(gsub("%", "", combined$mapping_efficiency))
combined$CHH_meth <- as.numeric(gsub("%", "", combined$CHH_meth))

#fig 4a
ggplot() +
  geom_violin(alpha = 0.45, fill = "#C7682D", data = combined, aes(species, mapping_efficiency, group = species)) +
  geom_boxplot(alpha = 0.65, fill = "#C7682D", width = 0.15, data = combined, aes(species, mapping_efficiency, group = species)) +
  geom_jitter(alpha = 0.75, color = "#C7682D", size = 6, data = combined, aes(species, mapping_efficiency, group = species)) +
  theme_classic(base_size = 40) +
  theme(axis.title.x = element_blank()) +
  ylab("Mapping efficiency (%)") +
  #labs(title = "Mapping Efficiency by NHP Species") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Capuchin", "Gelada", "Macaque")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(40,90)



## fig 4b
sites <- data.frame(observed = c(2349027, 3818703, 3546947),
                    expected = c(3343133, 5387280, 5486073),
                    species = c("capuchin", "gelada", "macaque"))
#barplot
sub_sites <- sites[,1:3]
sub_sites_long <- pivot_longer(sub_sites, cols = c("observed", "expected"), 
                               names_to = "variable", values_to = "value")

options(scipen = 999)
ggplot(sub_sites_long, aes(x = species, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic(base_size = 40) +
  labs(y = "Number of Sites") +
  scale_x_discrete(limits = c("macaque", "gelada", "capuchin"), labels = c("Macaque", "Gelada", "Capuchin")) +
  scale_fill_manual(values = c("expected" = "#C7682D", "observed" = "#FEB359")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")
