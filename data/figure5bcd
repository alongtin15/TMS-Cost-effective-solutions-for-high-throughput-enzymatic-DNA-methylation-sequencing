rrbs <- read.delim("model_summaries_RRBS.txt")
tms <- read.delim("model_summaries_TMS.txt")

tissues <- c("adrenal","heart","kidney","liver","lung","spleen")

rrbs$qval <- qvalue::qvalue(rrbs$pval)$qvalues
tms$qval <- qvalue::qvalue(tms$pval)$qvalues
#

## Assess correlation of effects --------
# overall correlation in effect sizes
cor.test(rrbs$beta, tms$beta)

# correlations per tissue - use this
correlation_df <- do.call(rbind, lapply(tissues, FUN = function(i) { 
  rrbs_tissue <- rrbs[which(rrbs$tissue==i),]
  tms_tissue <- tms[which(tms$tissue==i),]
  
  cor_test <- cor.test(rrbs_tissue$beta, tms_tissue$beta)
  data.frame(tissue = i,
             correlation = cor_test$estimate,
             pvalue = cor_test$p.value)
}))


## Overlap - fishers exact test --------
odds_ratio_df <- do.call(rbind, lapply(tissues, FUN = function(i) { 
  rrbs_tissue <- rrbs[which(rrbs$tissue==i),]
  tms_tissue <- tms[which(tms$tissue==i),]
  
  # Identify significant rows (pval < 0.05)
  rrbs_significant <- rrbs_tissue[which(rrbs_tissue$qval < 0.05),]$CpG
  tms_significant <- tms_tissue[which(tms_tissue$qval < 0.05),]$CpG
  
  # Determine overlaps and non-overlaps
  both_significant <- length(intersect(rrbs_significant, tms_significant))
  only_rrbs <- length(setdiff(rrbs_significant, tms_significant))
  only_tms <- length(setdiff(tms_significant, rrbs_significant))
  neither_significant <- nrow(rrbs_tissue) - (both_significant + only_rrbs + only_tms)
  
  # Fishers test 
  out_ft <- fisher.test(rbind(c(both_significant, only_rrbs), c(only_tms, neither_significant)))
  out_ft <- data.frame(pval = out_ft$p.value, 
                       estimate = out_ft$estimate, 
                       ci_1 = out_ft$conf.int[1], 
                       ci_2 = out_ft$conf.int[2], 
                       tissue = i,
                       n_both_significant = both_significant)
  return(out_ft)
}))


## Plots ---------
tmp<-dplyr::inner_join(rrbs[which(rrbs$tissue=="liver"),c("beta","CpG")] %>% 
        dplyr::rename("RRBS_beta"="beta"), 
  tms[which(tms$tissue=="liver"),c("beta","CpG")] %>% 
    dplyr::rename("TMS_beta"="beta"), by = "CpG")
op <- par(mar = c(5,7,4,2) + 0.1)
smoothScatter(x=tmp$RRBS_beta,y=tmp$TMS_beta,
              xlab = "RRBS effect size (liver-specific methylation)", 
              ylab= "TMS effect size (liver-specific methylation)",
              cex.lab=1.5, #change font size of axis labels
              cex.axis=1.5) #change font size of axis text   
par(op)
rm(tmp)
# save 5x5

tmp1 <- correlation_df %>% 
  mutate(tissue = stringr::str_to_title(tissue), 
         tissue = factor(tissue, levels = c("Adrenal","Liver","Spleen","Kidney","Heart","Lung"))) %>% 
  ggplot(aes(x = tissue, y = correlation)) + #, label = round(correlation, 2)
  geom_bar(stat = "identity", color = "black", fill = "grey80") +
  # geom_text(nudge_y = 0.05) +
  labs(x= "Tissue", y = "Correlation coefficent") + 
  scale_y_continuous(limits=c(0,1)) + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
rm(tmp1)

tmp2<-odds_ratio_df %>% 
  mutate(tissue = stringr::str_to_title(tissue), 
         tissue = factor(tissue, levels = rev(c("Adrenal","Liver","Spleen","Kidney","Heart","Lung")))) %>% 
  ggplot(aes(y = tissue, x = log2(estimate))) + #, label = n_both_significant
  geom_vline(xintercept = 0, lty=2, color = "grey50") + 
  geom_errorbar(aes(xmin = log2(ci_1), xmax = log2(ci_2))) + 
  geom_point() + 
  scale_x_continuous(limits = c(-1,4.5)) + 
  # geom_text(x = 4.5) + 
  labs(x = expression(Log[2]*" Fisher's Exact Test odds ratio", y = "")) + 
  theme_classic(base_size = 19) + theme(axis.title.y = element_blank())
