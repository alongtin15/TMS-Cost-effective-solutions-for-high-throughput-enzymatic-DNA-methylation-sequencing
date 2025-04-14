meanMeth_corr_data <- read.delim('meanMeth_corr_data.txt')

ggplot(meanMeth_corr_data %>% sample_n(30000, replace = FALSE),
       aes(x = meanMeth_rrbs, y = meanMeth_twist))+
  geom_point(col="#C7682D", alpha = 0.2, size = 4)+
  labs(x="RRBS", y="TTMS")+
  geom_abline(linetype="dashed", linewidth = 2)+
  theme_classic(base_size = 40)
