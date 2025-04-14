library(ggplot2)
library(tidyverse)

dragen <- read.delim('count_lines_dragen_22Aug2024.txt', header = F, sep = " ")
colnames(dragen) <- c("sampleID", "dragen_lines", "intersect", "intersect_unique", "twist_sub_dragen", "dragen_sub_twist")
dragen$sampleID <- gsub(".bam","",dragen$sampleID)
dragen$twist_probes <- 551803
dragen <- dragen[-grep("P1_NC",dragen$sampleID),]
dragen$experiment<-NA
for (i in c("NoME_65","65C_NoME","65C_2uLME","65C_4uLME","68C_NoME","68C_2uLME")) {
  dragen[grep(i,dragen$sampleID),]$experiment <- i
}
dragen[is.na(dragen$experiment),]$experiment <- "fragmentation"
dragen$experiment <- factor(dragen$experiment, levels = c("25_input","50_input","100_input","200_input","400_input",
                                                          "12plex", "24plex", "48plex","96plex",
                                                          "NoME_65","65C_NoME","65C_2uLME","65C_4uLME","68C_NoME","68C_2uLME",
                                                          'fragmentation'))
dragen$experiment <- gsub("NoME_65", "65C_NoME", dragen$experiment)
dragen <- dragen %>%
  filter(!experiment == "fragmentation")

## fig 2e
ggplot() + 
  geom_violin(data = dragen, alpha = 0.35, fill = "#71424E", aes(y = experiment, x = (dragen_sub_twist/dragen_lines)*100, group = experiment)) +
  geom_boxplot(data = dragen, alpha = 0.65, fill = "#71424E", width = 0.15, aes(y = experiment, x = (dragen_sub_twist/dragen_lines)*100, group = experiment)) + 
  geom_jitter(data = dragen, alpha = 0.75, color = "#71424E", size = 6, aes(y = experiment, x = (dragen_sub_twist/dragen_lines)*100, group = experiment)) + 
  theme_classic(base_size=40) +
  scale_y_discrete(limits = c("68C_2uLME","68C_NoME","65C_4uLME","65C_2uLME","65C_NoME"), 
                   labels = c("68C 2uL ME", "68C No ME", "65C 4uL ME", "65C 2uL ME", "65C No ME")) +
  xlab("% of Mapped Reads that are Off-Target") +
  #ylab("Experiment") +
  theme(axis.title.y = element_blank())


## fig 2f
ggplot() + 
  geom_violin(data = dragen, alpha = 0.35, fill = "#71424E", width = 1, aes(y = experiment, x = (intersect/twist_probes)*100, group = experiment)) +
  geom_boxplot(data = dragen, alpha = 0.65, fill = "#71424E", width = 0.15, aes(y = experiment, x = (intersect/twist_probes)*100, group = experiment)) + 
  geom_jitter(data = dragen, alpha = 0.65, color = "#71424E", size = 6, aes(y = experiment, x = (intersect/twist_probes)*100, group = experiment)) + 
  scale_y_discrete(limits = c("68C_2uLME","68C_NoME","65C_4uLME","65C_2uLME","65C_NoME"), 
                   labels = c("68C 2uL ME", "68C No ME", "65C 4uL ME", "65C 2uL ME", "65C No ME")) +
  xlab("% of Twist Probes Represented") +
  #ylab("Experiment") +
  xlim(65,100) +
  theme(axis.title.y = element_blank())
