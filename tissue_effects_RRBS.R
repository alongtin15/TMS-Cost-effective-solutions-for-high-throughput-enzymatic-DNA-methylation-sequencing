#!/usr/bin/env Rscript

### Model tissue effects in RRBS and TMS data for Amy's paper. Dec 30, 2024

library(dplyr)


# Read in data 
rrbs <- read.delim("/scratch/alongti2/EMseq/macaque/rrbs_common_pmeth.txt",header = T,sep = " ")
tms <- read.delim("/scratch/alongti2/EMseq/macaque/twist_common_pmeth.txt",header = T,sep = " ")
sample_to_tissue <- read.delim("/scratch/alongti2/EMseq/macaque/sample_to_tissue.csv",header = T,sep = ",")

# format 
sample_to_tissue$id <- gsub("-",".",gsub("11119-","",sample_to_tissue$Sample_Name))
sample_to_tissue <- sample_to_tissue %>% 
  mutate(adrenal = case_when(grantparent_tissueType == "adrenal" ~ 1, .default = 0), 
         heart = case_when(grantparent_tissueType == "heart" ~ 1, .default = 0),
         kidney = case_when(grantparent_tissueType == "kidney" ~ 1, .default = 0),
         liver = case_when(grantparent_tissueType == "liver" ~ 1, .default = 0),
         lung = case_when(grantparent_tissueType == "lung" ~ 1, .default = 0),
         spleen = case_when(grantparent_tissueType == "spleen" ~ 1, .default = 0),)

# tissues
tissues <- c("adrenal","heart","kidney","liver","lung","spleen")


## Remove hypo/hyper-methylated sites ---------------------
# hypo methylated sites are causing model failure, and removing nonvariable sites will also reduce testing burden 
rrbs_perc_meth <- data.frame(perc_meth = apply(rrbs, 1, FUN = median, na.rm=T), 
                             site = rownames(rrbs))
hist(rrbs_perc_meth$perc_meth,100)
rrbs_perc_meth_variable <- rrbs_perc_meth[which(rrbs_perc_meth$perc_meth > 0.1 & 
                                                  rrbs_perc_meth$perc_meth < 0.9),]

# variable sites in TMS dataset 
tms_perc_meth <- data.frame(perc_meth = apply(tms, 1, FUN = median, na.rm=T), 
                             site = rownames(tms))
hist(tms_perc_meth$perc_meth,100)
tms_perc_meth_variable <- rrbs_perc_meth[which(tms_perc_meth$perc_meth > 0.1 & 
                                                  tms_perc_meth$perc_meth < 0.9),]

## Get sites that are variably methylated for both technologies (they are highly correlated) 
dim(rrbs_perc_meth_variable)
dim(tms_perc_meth_variable)
length(which((rrbs_perc_meth_variable$site %in% tms_perc_meth_variable$site)==T))

variable_sites <- intersect(rrbs_perc_meth_variable$site, tms_perc_meth_variable$site)

# Filter for variable sites
rrbs_var <- rrbs[variable_sites,]
tms_var <- tms[variable_sites,]

rm(rrbs_perc_meth, tms_perc_meth, rrbs_perc_meth_variable, tms_perc_meth_variable, variable_sites)

#
## model - RRBS data ----------
rrbs_var_2 <- as.data.frame(t(rrbs_var)) # make rows = ids so you can merge with tissue type
rrbs_var_2$id <- rownames(rrbs_var_2)
rrbs_var_2 <- left_join(rrbs_var_2, sample_to_tissue[,c("id",tissues)], by = "id")

model_summaries_rrbs <- do.call(rbind, lapply(tissues, function(y) {
  do.call(rbind, lapply(names(rrbs_var_2)[-c(which(names(rrbs_var_2) %in% c("id", tissues)))], function(x) {
    rrbs_var_2_tmp <- rrbs_var_2[which(complete.cases(rrbs_var_2[[x]])), c(x, y)]
    colnames(rrbs_var_2_tmp) <- c("site","tissue")
    # 
    out <- summary(lm(site ~ tissue, data = rrbs_var_2_tmp))
    out <- as.data.frame(out$coefficients)[2,c(1,2,4)]
    colnames(out) <- c("beta","se","pval")
    out$CpG <- x
    out$tissue <- y
    out$n_samples <- nrow(rrbs_var_2_tmp)
    return(out)
  }))
}))

write.table(model_summaries_rrbs, "/scratch/mwatowic/tsimane_EMseq/tissue_effects_RRBSvTMS/model_summaries_RRBS.txt", quote = F, sep = "\t", row.names = F)

print("Successfully written all model outputs")

## module load r-4.2.2-gcc-11.2.0
## sbatch -p general -t 8:00:00 --mem=50GB tissue_effects_RRBS.R
