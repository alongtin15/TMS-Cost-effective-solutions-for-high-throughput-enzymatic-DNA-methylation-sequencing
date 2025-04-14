library(NatParksPalettes)
filt_twist_macaque_meth <- read.delim("filt_twist_macaque_meth.txt", header = T, sep = " ")

filt_twist_macaque_meth <- na.omit(filt_twist_macaque_meth)

#normalizing the data here
norm_meth <- scale(filt_twist_macaque_meth[,6:101])

#creating the correlation matrix
corr_matrix_meth <- cor(norm_meth)
ggcorrplot(corr_matrix_meth)

#running the PCA yay!
meth.pca <- princomp(corr_matrix_meth)
summary(meth.pca)

#making a scree plot to visualize the data
fviz_eig(meth.pca, addlabels = TRUE)

#extracting the loadings from the first three PCs
meth_loadings <- meth.pca$loadings[,1:3]
meth_loadings <- as.data.frame(meth_loadings)
setDT(meth_loadings, keep.rownames = "sampleID")

#using adding the tissue type to the loadings df
tissue <- macaque_meta[,4]
meth_loadings$tissue <- tissue

#plotting PC1 and PC2 coloring for tissue type
ggplot(meth_loadings, aes(x = Comp.1, y = Comp.2, color = tissue)) +
  geom_point(size = 6) +
  theme_classic(base_size = 40) +
  xlab("PC 1 (36.3%)") +
  ylab("PC 2 (23.5%)") +
  scale_color_manual(values = natparks.pals("Cuyahoga", n = 6))
