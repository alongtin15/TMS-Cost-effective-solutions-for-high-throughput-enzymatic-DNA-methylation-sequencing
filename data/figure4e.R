rrbs_common_pmeth <- read.delim('rrbs_common_pmeth.txt', sep = " ")
twist_common_pmeth <- read.delim('twist_common_pmeth.txt', sep = "")
rrbs_common_pmeth$avg_meth_rrbs <- apply(rrbs_common_pmeth, 1, function(x) mean(x, na.rm=TRUE))
twist_common_pmeth$avg_meth_twist <- apply(twist_common_pmeth, 1, function(x) mean(x, na.rm=TRUE))
combined <- merge(rrbs_common_pmeth, twist_common_pmeth, by = "row.names")

variable_sites <- combined %>%
  filter(avg_meth_rrbs > 0.1 & avg_meth_rrbs < 0.9 &
           avg_meth_twist > 0.1 & avg_meth_twist < 0.9)
dim(variable_sites)
  # 92692   195
rrbs_variable_sites <- variable_sites[,2:97]
twist_variable_sites <- variable_sites[,99:194]
colnames(rrbs_variable_sites) <- gsub(".x", "", colnames(rrbs_variable_sites))
colnames(twist_variable_sites) <- gsub(".y", "", colnames(twist_variable_sites))

rrbs_common_pmeth <- rrbs_common_pmeth[,1:96]
twist_common_pmeth <- twist_common_pmeth[,1:96]

results <- data.frame(sample = colnames(rrbs_common_pmeth))
results$beta5x_allsites <- NA
results$r.squared5x_allsites <- NA
results$beta5x_variable <- NA
results$r.squared5x_variable <- NA
for (i in 1:dim(results)[1]){
  rrbs_tmp <- as.data.table(rrbs_common_pmeth[, results$sample[i]])
  twist_tmp <- as.data.frame(twist_common_pmeth[, i])
  both1 <- cbind(rrbs_tmp, twist_tmp)
  results$beta5x_allsites[i]<-summary(lm(as.matrix(both1[,1]) ~ as.matrix(both1[,2])))$coefficients[2,1]
  results$r.squared5x_allsites[i]<-summary(lm(as.matrix(both1[,1]) ~ as.matrix(both1[,2])))$r.squared
  
  rrbs_tmp2 <- as.data.table(rrbs_variable_sites[, results$sample[i]])
  twist_tmp2 <- as.data.frame(twist_variable_sites[, i])
  both2 <- cbind(rrbs_tmp2, twist_tmp2)
  results$beta5x_variable[i] <- summary(lm(as.matrix(both2[,1]) ~ as.matrix(both2[,2])))$coefficients[2,1]
  results$r.squared5x_variable[i] <- summary(lm(as.matrix(both2[,1]) ~ as.matrix(both2[,2])))$r.squared
  print(i)
}

ggplot(results)+
  geom_density(position = "identity", alpha=0.75, aes(r.squared5x_allsites), fill = "#C7682D") +
  geom_density(position = "identity", alpha=0.75, aes(r.squared5x_variable), fill = "#FEB359") +
  labs(x="R-squared", y="Density")+
  theme_classic(base_size=40)+
  theme(legend.position = "none")
