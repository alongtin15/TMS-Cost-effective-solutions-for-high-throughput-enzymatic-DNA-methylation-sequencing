## Epigenetic age predictions for EPIC and TMS data - Amy paper. Jan 2024 
library(dplyr)
library(ggplot2)

### Methylation percentages from paired EPIC/TMS data ---------
# TMP <- readRDS("path/to/unique_EPIC_EMseq_combined_02jan25.rds")
# TMP$chr_loc <- paste0("chr",TMP$chromosome,"_",TMP$start)
# 
# ## Connect with Illumina probes 
# sites<-as.data.frame(TMP$chr_loc)
# colnames(sites) <- "chr_loc"
# sites$id<-1:dim(sites)[1]
# 
# ## Get EPIC probes and subset to just these sites; rename rows to probe names 
# probes=read.csv('path/to/Illumina_EPIC_manifests/infinium-methylationepic-v-1-0-b5-manifest-file_v2.csv')
# probes$chr_loc<-paste(probes$CHR_hg38,probes$Start_hg38+1,sep='_')
# both<-merge(sites,probes[,c('chr_loc','Name')],by='chr_loc')
# both<-both[order(both$id),]
# 
# epic<-TMP[both$id, c("chr_loc", colnames(TMP)[grep("App574", colnames(TMP))])]
# tms<-TMP[both$id, c("chr_loc", colnames(TMP)[grep("10274-AJ", colnames(TMP))])]
# sites_2<-sites[both$id,]
# 
# # Merge with probe names
# epic_2 <- merge(epic,probes[,c("chr_loc","IlmnID")], by = "chr_loc")
# rownames(epic_2) <- epic_2$IlmnID
# epic_2$chr_loc <- NULL
# epic_2$IlmnID <- NULL
# 
# tms_2 <- merge(tms,probes[,c("chr_loc","IlmnID")], by = "chr_loc")
# rownames(tms_2) <- tms_2$IlmnID
# tms_2$chr_loc <- NULL
# tms_2$IlmnID <- NULL
# # 
# rm(TMP, probes, both, sites, sites_2)
# 
# 
### Use PC clock package for epigenetic age estimate --------------
# clocksDir <- "path/to/PC-Clocks-main/" 
# source(paste(clocksDir, "run_calcPCClocks.R", sep = ""))
# source(paste(clocksDir, "run_calcPCClocks_Accel.R", sep = ""))
# 
# # Format data 
# epic_pcClocks <- t(epic_2)
# tms_pcClocks <- t(tms_2)
# 
# pheno_epic <- data.frame(SampleID = rownames(epic_pcClocks), 
#                          Female = 1,
#                          Age = 40)
# 
# pheno_tms <- data.frame(SampleID = rownames(tms_pcClocks),
#                         Female = 1,
#                         Age = 40)
# 
# # Get the PC Clocks values and the PC Clock Acceleration values
# PCClock_DNAmAge_epic <- calcPCClocks(path_to_PCClocks_directory = clocksDir,
#                                      datMeth = epic_pcClocks,
#                                      datPheno = pheno_epic)
# 
# PCClock_DNAmAge_tms <- calcPCClocks(path_to_PCClocks_directory = clocksDir,
#                                     datMeth = tms_pcClocks,
#                                     datPheno = pheno_tms)
# 
# # Remove grimAge and components b/c they use Age/Sex info and ours is dummy coded
# PCClock_DNAmAge_epic <- PCClock_DNAmAge_epic[,c("SampleID","PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCDNAmTL")]
# PCClock_DNAmAge_tms <- PCClock_DNAmAge_tms[,c("SampleID","PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCDNAmTL")]
# 
# # Rename Cai epic samples with same sample ID as TMS 
# cai_names <- read.delim("~/Downloads/cai_file_names.txt",header = F,sep = "\t")
# colnames(cai_names) <- c("tms_name","epic_name")
# PCClock_DNAmAge_epic <- merge(PCClock_DNAmAge_epic, cai_names, by.x = "SampleID", by.y = "epic_name")
# PCClock_DNAmAge_epic <- PCClock_DNAmAge_epic[,-1]
# colnames(PCClock_DNAmAge_epic)[colnames(PCClock_DNAmAge_epic) == "tms_name"] <- "SampleID"
# 
# # Combine clock estimates from two technologies 
# PCClock_DNAmAge_epic <- PCClock_DNAmAge_epic %>%
#   rename_with(~ paste0(., "_epic"), -SampleID)
# 
# PCClock_DNAmAge_tms <- PCClock_DNAmAge_tms %>%
#   rename_with(~ paste0(., "_tms"), -SampleID)
# 
# # Combine the data frames by id
# pcClocks_comb <- full_join(PCClock_DNAmAge_tms, PCClock_DNAmAge_epic, by = "SampleID")
# 
# # Write file 
# write.table(pcClocks_comb, "path/to/pcClock_EPIC_TMS_6Jan2024.txt", quote = F, sep = "\t", row.names = F)


### Compare clock estimated ages and make figure -------------
pcClocks_comb <- read.delim("~/Desktop/pcClock_Cai_EPIC_TMS_6Jan2024.txt",header = T,sep = "\t")

cor_results <- do.call(rbind, lapply(c("PCHorvath1","PCHorvath2","PCHannum","PCPhenoAge","PCDNAmTL"), function(i) {
  cor_test <- cor.test(pcClocks_comb[[paste0(i,"_epic")]], pcClocks_comb[[paste0(i,"_tms")]])
  data.frame(Clock = i,
             Correlation = cor_test$estimate,
             Pvalue = cor_test$p.value)
}))
cor_results

tmp<-cor_results %>%
  mutate(Clock = factor(Clock, levels = c(cor_results$Clock))) %>% 
  ggplot(aes(x = Clock, y = Correlation)) + #, label = round(Correlation, 2)
  geom_bar(stat = "identity", color = "black", fill = "grey80") +
  # geom_text(nudge_y = 0.05) + 
  labs(x= "Clock", y = "Correlation coefficent") + 
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(tmp, file = "~/Documents/_Vanderbilt/Tsimane/DNA_methylation/EMseq_pilot_2023/epigenetic_clock_age_predictions/clock_corr.pdf", device = "pdf",width = 5, height = 5)


# Plot of Hannum estimates (blood-based clock)
pcClocks_comb %>% 
  ggplot(aes(x = PCHannum_epic, y = PCHannum_tms)) + 
  geom_point() + 
  theme_classic()

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
  theme_classic() + 
  facet_wrap(.~clock, scales="free", nrow=1)
