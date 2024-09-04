##### This script compares the results of TWIST-CAPTURE and RRBS procedures.

##### Copy and modified 29/07/2024 by Baptiste Sadoughi

rm(list = ls())

setwd("/path/to/working/directory/")

library_list <- c("tidyverse","comethyl","GenomicRanges")
lapply(library_list, require, character.only=TRUE)

# at the moment this is causing massive issues and useNames.NA should be set to warning instead of error
options(matrixStats.useNames.NA = "deprecated")


#------------------------------------------------------------------------------
########################################### PREPARE DATA
#------------------------------------------------------------------------------

# Load Twist capture rhesus macaque
twist <- readRDS("/path/to/twist/RDS/file/")

# Simplify sample names
colnames(twist) <- gsub(pattern = "XX","",
                        str_split_i(colnames(twist),"/",7))
colnames(twist) <- gsub(pattern = "XX","",colnames(twist))

# Load metadata
metadata_lid = read.table("/path/to/metadata/table/", sep = "\t", header = TRUE)

# Load Twist sample DIDs
identifiers = readxl::read_xlsx("/path/to/files/matching/sampleIDs/to/fileIDs/", col_names = TRUE)

# clean headers
colnames(identifiers) <- gsub(" ","_",colnames(identifiers))

# Merge metadata with Twist data
metadata <- identifiers %>%
  left_join(metadata_lid %>%
              select(parentName, monkey_id,individual_sex, age_at_sampling, grantparent_tissueType),
            join_by("Description" == "parentName")) %>%
  unique() %>% mutate(Sample_Name = gsub("XX","",Sample_Name)) %>% filter(complete.cases(.))

# Load rrbs
rrbs=readRDS("/path/to/RRBS/RDS/file/")

# Simplify sample names
colnames(rrbs) <- gsub("XX", "",
                       str_split_i(colnames(rrbs), "/", 6))

# Subset rrbs to match Twist
for_comparison <- metadata$Description
rrbs_metadata <- metadata_lid %>% filter(<sample name df1> %in% <sample name df2>)

# Keep only one matching rrbs
rrbs_metadata <- rrbs_metadata %>%
  group_by(<sample ID>, <tissue type>) %>%
  mutate(n = n_distinct(<sample identifier>), max_unique = max(unique))#keep the sample with the highest unique reads

rrbs_metadata <- rrbs_metadata %>%
  filter(n == 1 | (n > 1 & unique == max_unique)) %>%
  dplyr::select(-n, -max_unique) %>%
  ungroup()

rrbs_metadata <- rrbs_metadata %>% select(<sample identifier>)
metadata <- merge(metadata,rrbs_metadata, by.x = "<sample identifier df1>", by.y = "<sample identifier df2>")

rrbs <- rrbs[,metadata$<sampleID>]
twist <- twist[,metadata$<sampleID>]

rm(metadata_lid,rrbs_metadata)


#------------------------------------------------------------------------------
########################## FILTERING FOR COVERAGE
#------------------------------------------------------------------------------

# How many sites are covered with each technique?

# Remove site with zero coverage
num_cores = 21
chrs=paste0("",c(1:20,"X")) #adapt to genome of the species

rrbs_cov1_list <- parallel::mclapply(chrs,function(x){
  #select one chr at a time
  chr <- chrSelectBSseq(rrbs, seqnames = x, order = TRUE)
  
  # apply filtering
  rrbs_chr <- filterCpGs(chr,
                         cov = 5, perSample = 0.75, verbose = FALSE,
                         save = FALSE, file = NULL)
  return(rrbs_chr)
}, mc.cores = num_cores)

twist_cov1_list <- parallel::mclapply(chrs,function(x){
  #select one chr at a time
  chr <- chrSelectBSseq(twist, seqnames = x, order = TRUE)
  
  # apply filtering
  twist_chr <- filterCpGs(chr,
                          cov = 5, perSample = 0.75,verbose = FALSE,
                          save = FALSE, file = NULL)
  return(twist_chr)
}, mc.cores = num_cores)

# Add a percent methylation matrix on CpG sites
rrbs_cov1_list=parallel::mclapply(
  rrbs_cov1_list,function(x){
    pmeth=getCoverage(x,type="M")/getCoverage(x,type="Cov")
    obj=x
    obj@assays@data@listData[["pmeth"]]=pmeth
    return(obj)
  },mc.cores = num_cores)

twist_cov1_list=parallel::mclapply(
  twist_cov1_list,function(x){
    pmeth=getCoverage(x,type="M")/getCoverage(x,type="Cov")
    obj=x
    obj@assays@data@listData[["pmeth"]]=pmeth
    return(obj)
  },mc.cores = num_cores)

# Mean percent methylation per row
chrs=c(1:21)
meanMeth_rrbs_list <- parallel::mclapply(chrs,function(x){
  meanMeth <- rowMeans(rrbs_cov1_list[[x]]@assays@data@listData[["pmeth"]], na.rm = TRUE)
  
  return(meanMeth)
}, mc.cores = (length(chrs)+1))

meanMeth_rrbs <- unlist(meanMeth_rrbs_list)

meanMeth_twist_list <- parallel::mclapply(chrs,function(x){
  meanMeth <- rowMeans(twist_cov1_list[[x]]@assays@data@listData[["pmeth"]], na.rm = TRUE)
  
  return(meanMeth)
}, mc.cores = (length(chrs)+1))

meanMeth_twist <- unlist(meanMeth_twist_list)

# Create a data frame for the two vectors
df <- data.frame(
  value = c(meanMeth_rrbs, meanMeth_twist),
  group = factor(rep(c("meanMeth_rrbs", "meanMeth_twist"), times = c(length(meanMeth_rrbs),length(meanMeth_twist))))
)


#------------------------------------------------------------------------------
##################################### INTERSECTING THE BSSEQ
#------------------------------------------------------------------------------

###### GRangesList
######
# FUNCTION to convert bsseq list to GRangesList
bsseq_to_grl <- function(bsseq_list) {
  gr_list <- lapply(bsseq_list, function(bsseq) {
    chrom <- seqnames(bsseq)
    start <- start(bsseq)
    end <- end(bsseq)
    gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
    return(gr)
  })
  return(GRangesList(gr_list))
}

extract_common_cpgs <- function(bsseq_list_1, bsseq_list_2,bsseq_common_1_name = "bsseq_list_1_common", bsseq_common_2_name = "bsseq_list_2_common") {
  # Create GRangesList from bsseq_list_1
  bsseq_grl_1 <- bsseq_to_grl(bsseq_list_1)
  
  # Create GRangesList from bsseq_list_2
  bsseq_grl_2 <- bsseq_to_grl(bsseq_list_2)
  
  # Create an empty GRangesList to store common CpGs
  common_cpg <- GRangesList()
  
  # Iterate through each chromosome and find the common CpG sites
  for (chr in seq_along(bsseq_list_1)) {
    # Find common CpG sites between bsseq_grl_1 and bsseq_grl_2 for the current chromosome
    common_sites <- GenomicRanges::intersect(bsseq_grl_1[[chr]], bsseq_grl_2[[chr]])
    
    # Create a logical vector indicating which elements in bsseq_grl_1[[chr]] overlap with common_sites
    overlap_logical <- findOverlaps(bsseq_grl_1[[chr]], common_sites)
    # Extract the indices of overlapping ranges in bsseq_grl_1[[chr]]
    overlap_indices <- queryHits(overlap_logical)
    # Subset bsseq_grl_1[[chr]] based on the overlap_indices
    subset_bsseq_1 <- bsseq_grl_1[[chr]][overlap_indices]
    
    # Assign Subsetted Ranges to common_cpg
    common_cpg[[chr]] <- subset_bsseq_1
  }
  
  # Finally subset the bsseq objects to keep common CpGs
  # Initialize lists to store the bsseq objects
  bsseq_common_1 <- list()
  bsseq_common_2 <- list()
  
  # Iterate through each chromosome in common_cpg
  for(chr in seq_along(bsseq_list_1)) {
    # Get the common CpG sites
    common_cpg_chr <- common_cpg[[chr]]
    
    # Subset the bsseq lists based on common CpG sites
    bsseq_common_1[[chr]] <- subsetByOverlaps(bsseq_list_1[[chr]], common_cpg_chr)
    bsseq_common_2[[chr]] <- subsetByOverlaps(bsseq_list_2[[chr]], common_cpg_chr)
  }
  
  # Save the output lists as separate objects within the environment
  assign(bsseq_common_1_name, bsseq_common_1, envir = .GlobalEnv)
  assign(bsseq_common_2_name, bsseq_common_2, envir = .GlobalEnv)
}

# Apply function to the bsseq objects
extract_common_cpgs(bsseq_list_1=rrbs_cov1_list, bsseq_list_2=twist_cov1_list, bsseq_common_1_name = "rrbs_common", bsseq_common_2_name = "twist_common")

#-------------------------------------------------------------------------------
############################### TWIST RRBS COMPARISONS
#-------------------------------------------------------------------------------

##### COMPARE DISTRIBUTION OF PERCENT METHYLATION

# Mean percent methylation per row
chrs=c(1:21)
meanMeth_rrbs_list <- parallel::mclapply(chrs,function(x){
  meanMeth <- rowMeans(rrbs_common[[x]]@assays@data@listData[["pmeth"]], na.rm = TRUE)
  
  return(meanMeth)
}, mc.cores = (length(chrs)+1))

meanMeth_rrbs <- unlist(meanMeth_rrbs_list)

meanMeth_twist_list <- parallel::mclapply(chrs,function(x){
  meanMeth <- rowMeans(twist_common[[x]]@assays@data@listData[["pmeth"]], na.rm = TRUE)
  
  return(meanMeth)
}, mc.cores = (length(chrs)+1))

meanMeth_twist <- unlist(meanMeth_twist_list)

# Create a data frame for the two vectors
df <- data.frame(
  value = c(meanMeth_rrbs, meanMeth_twist),
  group = factor(rep(c("meanMeth_rrbs", "meanMeth_twist"), each = length(meanMeth_rrbs)))
)

# Correlation of meanMeth for paired CpG
meanMeth_corr_data <- as.data.frame(cbind(meanMeth_rrbs,meanMeth_twist))
cor.meanMeth.1<-with(meanMeth_corr_data,cor.test(meanMeth_rrbs,meanMeth_twist, method="pearson"))

random_rrbs <- sample(meanMeth_corr_data$meanMeth_rrbs)
meanMeth_corr_data <- as.data.frame(cbind(meanMeth_corr_data,random_rrbs))
cor.meanMeth.1.rand<-with(meanMeth_corr_data,cor.test(meanMeth_twist,random_rrbs, method="pearson"))

#------------------------------------------------------------------------------
################## FUNCTIONS TO COMPARE MATCHED CPG AND MATCHED SAMPLES
#------------------------------------------------------------------------------
# below I wrote a number of functions to compare percent methylation across matched CpG sites,
# or across matched samples.
# Note : this is ran on datasets for which only CpG sites in common are left. But it is still necessary to
# restrict to CpG sites which are covered for a given matched pair of samples.
# Note also that row-wise comparison can be performed on a list, whereas for column-wise (i.e. sample) a whole bsseq is needed.

# Function to combine data from all chromosomes for each sample
combine_bsseq <- function(bsseq_list) {
  combined_data <- do.call(rbind, lapply(bsseq_list, function(bsseq) {
    bsseq@assays@data@listData[["pmeth"]]  # Assuming this extracts the relevant matrix
  }))
  return(combined_data)
}

combined_rrbs <- combine_bsseq(rrbs_common)  # Create combined data for RRBS
combined_twist <- combine_bsseq(twist_common)  # Create combined data for TWIST

# Calculate means per row for combined data (this takes time!)
means_rrbs <- rowMeans(combined_rrbs, na.rm = TRUE)
means_twist <- rowMeans(combined_twist, na.rm = TRUE)

num_samples <- ncol(combined_rrbs)  # Number of samples based on columns of the combined matrix

# Function to compute linear regression R² for combined data and variable sites
compute_combinedbsseq_regressionR2 <- function(combined_bsseq1, combined_bsseq2, means_bsseq1, means_bsseq2) {
  
  combined_bsseq1 <- as.matrix(combined_bsseq1)
  combined_bsseq2 <- as.matrix(combined_bsseq2)
  
  # Ensure there are no NA values
  non_na_indices <- complete.cases(combined_bsseq1, combined_bsseq2)
  
  non_na_bsseq1 <- combined_bsseq1[non_na_indices, , drop = FALSE]
  non_na_bsseq2 <- combined_bsseq2[non_na_indices, , drop = FALSE]
  
  # Identify variable sites based on the pre-computed means
  variable_indices <- (means_bsseq1 > 0.1 & means_bsseq1 < 0.9) & (means_bsseq2 > 0.1 & means_bsseq2 < 0.9)
  
  non_variable_indices <- !variable_indices
  
  # Subset variable_indices based on non-NA indices
  variable_indices <- variable_indices[non_na_indices]
  non_variable_indices <- non_variable_indices[non_na_indices]
  
  # Calculate the non-variable R²
  if (any(non_variable_indices) && sum(non_variable_indices) >= 2) {
    r.squared.nonvar <- summary(lm(non_na_bsseq1[non_variable_indices, , drop = FALSE] ~ non_na_bsseq2[non_variable_indices, , drop = FALSE]))$r.squared
  } else {
    r.squared.nonvar <- NA
  }
  
  # Calculate R² for variable sites
  if (any(variable_indices) && sum(variable_indices) >= 2) {
    r.squared.var <- summary(lm(non_na_bsseq1[variable_indices, , drop = FALSE] ~ non_na_bsseq2[variable_indices, , drop = FALSE]))$r.squared
  } else {
    r.squared.var <- NA
  }
  
  return(list(R_squared_nonvar = r.squared.nonvar, R_squared_var = r.squared.var))
}

results_linearR2 <- list()  # Initialize an empty list to store results

results_linearR2 <- mclapply(1:num_samples, function(idx) {
  compute_combinedbsseq_regressionR2(
    combined_rrbs[, idx, drop = FALSE],  # Access the specific column for the sample
    combined_twist[, idx, drop = FALSE],
    means_rrbs,
    means_twist
  )
}, mc.cores = 30)  # Adjust the number of cores as needed

# Print the results
print(results_linearR2)

# Convert the list of results into a data frame with Sample_Index included
regressionR2_coefficients <- do.call(rbind, lapply(seq_along(results_linearR2), function(idx) {
  res <- results_linearR2[[idx]]
  data.frame(Sample_Index = idx, R_squared_nonvar = res$R_squared_nonvar, R_squared_var = res$R_squared_var)
}))

###################### END