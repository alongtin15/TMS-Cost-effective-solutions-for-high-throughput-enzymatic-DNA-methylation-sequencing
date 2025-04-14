## loading libraries and data ----
library(tidyverse)
library(data.table)
library(ggplot2)
library(corrr)

pcClocks_comb <- read.delim("pcClock_Cai_EPIC_TMS_6Jan2024.txt",header = T,sep = "\t")

cor_results <- do.call(rbind, lapply(c("PCHorvath1","PCHorvath2","PCHannum","PCPhenoAge","PCDNAmTL"), function(i) {
  cor_test <- cor.test(pcClocks_comb[[paste0(i,"_epic")]], pcClocks_comb[[paste0(i,"_tms")]])
  data.frame(Clock = i,
             Correlation = cor_test$estimate,
             Pvalue = cor_test$p.value)
}))
cor_results

## plot for main text
cor_results %>%
  mutate(Clock = factor(Clock, levels = c(cor_results$Clock))) %>% 
  ggplot(aes(x = Clock, y = Correlation, label = round(Correlation, 2))) +
  geom_bar(stat = "identity", color = "white", fill = "#396EC2") +
  geom_text(nudge_y = 0.05) + 
  labs(x= "Clock", y = "Correlation coefficent") + 
  theme_classic(base_size = 40) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())

# Plot of Hannum estimates (blood-based clock)
pcClocks_comb %>% 
  ggplot(aes(x = PCHannum_epic, y = PCHannum_tms)) + 
  geom_point(size = 5, alpha = 0.75) + 
  theme_bw(base_size = 40) +
  xlab("EPIC") +
  ylab("TMS") +
  labs(title = "Hannum Estimates")
summary(lm(pcClocks_comb$PCHannum_tms ~ pcClocks_comb$PCHannum_epic))
  #R-squared: 0.9118

# Plot of all clock comparisons
pcClocks_comb %>% 
  pivot_longer(cols = c("PCHorvath1_epic","PCHorvath2_epic","PCHannum_epic","PCPhenoAge_epic","PCDNAmTL_epic",
                        c("PCHorvath1_tms","PCHorvath2_tms","PCHannum_tms","PCPhenoAge_tms","PCDNAmTL_tms")), 
               values_to = "clock_estimate", names_to = "clock_type") %>% 
  mutate(clock = gsub("\\_.*", "", clock_type),
         type = gsub(".*_", "", clock_type)) %>% 
  dplyr::select(c("SampleID","clock","type","clock_estimate")) %>% 
  pivot_wider(names_from = "type",values_from = "clock_estimate") %>% 
  mutate(clock = factor(clock, levels = c("PCHorvath1","PCHorvath2","PCHannum","PCPhenoAge","PCDNAmTL"))) %>% 
  ggplot(aes(x = epic, y = tms)) + 
  geom_point() + 
  theme_classic(base_size = 20) + 
  xlab("EPIC") +
  ylab("TMS") +
  facet_wrap(.~clock, scales="free", nrow=1)
