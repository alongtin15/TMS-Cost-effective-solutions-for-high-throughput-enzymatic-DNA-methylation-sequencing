library(ggplot2)
combined_avg_meth <- read.delim('combined_avg_meth.txt')

## plotting ----
ggplot(combined_avg_meth, aes(avg_meth*100, species, group = species)) +
  geom_density_ridges(fill = "#C7682D", alpha = 0.75) +
  theme_classic(base_size = 40) +
  xlab("Average Methylation (%)") +
  ylab("Species") +
  scale_y_discrete(limits = c("capuchin", "gelada", "macaque"),
                   labels = c("Capuchin", "Gelada", "Macaque")) +
  theme(legend.position = "none")
