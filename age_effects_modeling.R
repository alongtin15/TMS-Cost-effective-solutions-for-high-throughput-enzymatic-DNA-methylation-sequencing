### Beta Binomial models of age effects
library(aod)
library(dplyr)
### Read in data -----------------------------
## Counts and methylated counts 
total_counts <- read.delim("total_counts.txt", header = T, sep = "\t",check.names = F)
mcounts <- read.delim("mcounts.txt", header = T, sep = "\t",check.names = F)
rownames(total_counts)<-total_counts$site
rownames(mcounts)<-mcounts$site
total_counts$site<-NULL
mcounts$site<-NULL
## Metadata 
meta_EMseq_pilot <- read.delim("meta.csv", header = T,sep = ",")
## Cell counts 
cell_counts <- read.delim("cell_counts.txt",header = T,sep = "\t")
cell_counts <- cell_counts[,c("B","CD4T","CD8T")]
cell_counts <- cell_counts[-95,]
# Bind metadata and cell counts 
meta_EMseq_pilot <- cbind(meta_EMseq_pilot, cell_counts)
### Model -----------------------------
print("starting modeling")
out_96 <- do.call(rbind, lapply(row_start:row_end, function(x) { 
  
  print(paste0("starting model: ", x))
  
  # Merge counts and metadata - for this modeling fxn, counts need to be a column in the metadata 
  counts_df <- merge(as.data.frame(t(mcounts_new[x,])) %>% tibble::rownames_to_column(var = "id"), 
                     as.data.frame(t(total_counts[x,])) %>% tibble::rownames_to_column(var = "id"), 
                     by = c("id"))
  colnames(counts_df) <- c("id","mcounts","total_counts")
  counts_df <- inner_join(counts_df, 
                          meta_EMseq_pilot[,c("sample_id_short", "male", "age","B","CD4T","CD8T")],
                          by=c("id" = "sample_id_short")) 
 
  # Remove NAs
  counts_df <- counts_df[complete.cases(counts_df$male) & complete.cases(counts_df$age),]
  # Model
  mod_out<-summary(aod::betabin(cbind(mcounts, total_counts - mcounts) ~ 
                                  scale(age) + scale(male) + scale(B) + scale(CD4T) + scale(CD8T), ~ 1, data = counts_df))
  
    # Format output for each variable of interest
    age_out <- mod_out@Coef[2,c(1,4)]
    colnames(age_out) <- paste0(colnames(age_out)," ",rownames(age_out))
    sex_out <- mod_out@Coef[3,c(1,4)]
    colnames(sex_out) <- paste0(colnames(sex_out)," ",rownames(sex_out))
    B_out <- mod_out@Coef[4,c(1,4)]
    colnames(B_out) <- paste0(colnames(B_out)," ",rownames(B_out))
    CD4T_out <- mod_out@Coef[5,c(1,4)]
    colnames(CD4T_out) <- paste0(colnames(CD4T_out)," ",rownames(CD4T_out))
    CD8T_out <- mod_out@Coef[6,c(1,4)]
    colnames(CD8T_out) <- paste0(colnames(CD8T_out)," ",rownames(CD8T_out))
    
    # Combine model outputs 
    out <- cbind(age_out, sex_out, B_out, CD4T_out, CD8T_out)
    out$n <- nrow(counts_df)
    out$site <- rownames(mcounts[x,])
    return(out)
}))
colnames(out_96)
colnames(out_96)[1:10] <- c("beta_age","pval_age","beta_sex","pval_sex",
                           "beta_B","pval_B","beta_CD4T","pval_CD4T","beta_CD8T","pval_CD8T")
### save output
write.csv(out_96, paste0("out_",Sys.Date(),"_", formatC(thischunk, width = 6, flag = 0), ".csv"))
### Print info 
print(paste0("n(models) = ", nmodels,
             "\nn(chunks) = ", nchunks,
             "\nchunk #", thischunk,
             "\nchunk size = ", chunksize,
             "\n start at row: ", row_start,
             "\n end at row: ", row_end))
quit(save="no")