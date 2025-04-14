library(ggplot2)
library(tidyverse)

median_prop_tsi_l <- read.delim('median_prop_tsi_long.txt')
median_prop_tsi_l %>% 
  filter(yes_no=="1") %>% 
  ggplot(aes(y = reorder(state_name,median_proportion), x = (median_proportion)*100)) + 
  ggridges::geom_density_ridges(show.legend = F, fill="#123067", alpha = 0.5) +
  labs(x="Median DNA Methylation Level (%)", y="") + 
  theme_classic(base_size = 40) + 
  theme(axis.title.y = element_blank()) + 
  scale_x_continuous(limits=c(0,100))
  #save as 12x15 landscape

cpgs_per_region_tsi_96plex <- read.delim('cpgs_per_region_tsi_96plex.txt')
state_order <- c("Quiescent/Low","Strong transcription","Heterochromatin","Weak Repressed PolyComb",
                 "Weak transcription","Genic enhancers","Gene body",
                 "Transcr. at gene 5' and 3'","Promoters","ZNF genes & repeats","Enhancers",
                 "Repressed PolyComb","Flanking Active TSS","Bivalent Enhancer","Active TSS",
                 "Bivalent/Poised TSS","Flanking Bivalent TSS/Enh")
cpgs_per_region_tsi_96plex %>%
  mutate(state_name = factor(state_name, levels = rev(state_order))) %>%
  ggplot(aes(x = state_name)) + 
    geom_bar(fill = "#123067", width = 0.85) + 
    labs(x="",y="Count") + 
    coord_flip() + 
    theme_classic(base_size = 40) + 
    scale_y_continuous(labels = scales::comma)
