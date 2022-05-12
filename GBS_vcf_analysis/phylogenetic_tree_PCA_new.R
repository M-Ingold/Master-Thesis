library(ape)
library(vcfR)
library(cluster)
library(ggplot2)
library(factoextra)
library(dendextend)
library(readxl)
library(dplyr)
library(ggfortify)
library(ggrepel)
library(ggpubr)
library(poppr)

vcfTetra <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf")
vcfDiploid <- read.vcfR("../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf")

dosagetetra <- extract.gt(vcfTetra, as.numeric = F)

# could be done in fewer lines, but works for now
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/0", 0)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/1", 1)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/1/1", 2)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/1/1/1", 3)
dosagetetra <- replace(dosagetetra, dosagetetra == "1/1/1/1", 4)
# to-do: impute NAs
#dosagetetra[is.na(dosagetetra)] <- 0
dosagetetra <- as.numeric(na.omit(dosagetetra))

dosagetetra <- as.numeric(dosagetetra)

matrix_tetra <- t(matrix(unlist(dosagetetra), nrow = length(dosagetetra)/261, ncol = 261))


# get diploid dosages. Alternatively one could diploidize matrix_tetra
dosagedip <- extract.gt(vcfDiploid, as.numeric = F)

dosagedip <- replace(dosagedip, dosagedip == "0/0", 0)
dosagedip <- replace(dosagedip, dosagedip == "0/1", 1)
dosagedip <- replace(dosagedip, dosagedip == "1/1", 2)
#dosagedip[is.na(dosagedip)] <- 0
dosagedip <- as.numeric(na.omit(dosagedip))

matrix_dip <- t(matrix(unlist(dosagedip), nrow = length(dosagedip)/261, ncol = 261))




# import potato genotype data from excel sheet
GBS02 <- read_excel("../data/GBS2/GBS02_sample-IDs_14.04.21.xlsx", skip = 6)
GBS02 <- rename(GBS02, Identifier = "Sample identifier")
GBS01 <- read_excel("../data/GBS1/22-01-2021_Zuordnung_genotypes_SB.aktualisiert.xlsx")

# poorly merge datasets, works for now
intersect <- intersect(names(GBS01), names(GBS02))
GBS <- merge(GBS01, GBS02, by=intersect, all = T)

df <- data.frame(Identifier=gsub("Sample_", "",colnames(vcfTetra@gt)[-1]))
Samples <- merge(df, GBS, by="Identifier")

# add proper breed names to rows
rownames(matrix_tetra) <- Samples$VARIETY
rownames(matrix_dip) <- Samples$VARIETY


# # function to create tanglegrams. Gives an error message about labels not matching, no alternative to copy-pasting :(
# tangle <- function(dend1,dend2){
#   dend_list <- dendlist(dend1,dend2)
#   tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
#              common_subtrees_color_lines = T, # Turn-off line colors 
#              common_subtrees_color_branches = TRUE, # Color common branches 
#              main = paste("entanglement =", round(entanglement(dend_list, 2)))
#   )
# }

# load ADMIXTURE data
admixture=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.5.Q")
rownames(admixture) <- Samples$VARIETY

# assign population membership by highest membership percentage
admixture$subpopulation <- max.col(admixture)

# mark breeds with low membership as admixed
# optional
for (i in 1:length(admixture$V1)){
  if (max(admixture[i,1:5])<0.8) {
    admixture[i,6] <- 6
  }
}

# to avoid having to copy-paste a lot of plots, the following is optional. Just comment it out if subsetting by subpopulation membership isn't needed
#########################################################
min <- 0.8
admixture <- subset(admixture, V1 > min|V2 > min|V3 > min|V4 > min|V5 > min)

# subset for non-admixed breeds
matrix_tetra <- matrix_tetra[rownames(admixture), ]
matrix_dip <- matrix_dip[rownames(admixture), ]

# after subsetting, some rows all contain the same number. These have to be removed for PCA, in this case using the variance.
matrix_tetra <- (matrix_tetra)[ , which(apply((matrix_tetra), 2, var) != 0)]
matrix_dip <- (matrix_dip)[ , which(apply((matrix_dip), 2, var) != 0)]

##########################################################
admixture[,6] <- as.character(admixture[,6])


# calculate distances between samples
nei_dist_dip <- nei.dist((matrix_dip))
nei_dist_tetra <- nei.dist((matrix_tetra))

nei_dist_dip_wardd2 <- (hclust(d= nei_dist_dip, method = "ward.D2"))
nei_dist_tetra_wardd2 <- (hclust(d= nei_dist_tetra, method = "ward.D2"))

# create tanglegram comparing diploid and tetraploid data
dend1 <- as.dendrogram(nei_dist_dip_wardd2) 
dend2 <- as.dendrogram(nei_dist_tetra_wardd2)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = "Diploid vs. tetraploid, WardD2 using Nei's distance",
           cex_main = 1.5,
           sub = paste("entanglement =", round(entanglement(dend_list), 2)),
           cex_sub = 1,
           lab.cex = 0.5
)


# calculate distances between samples
euc_dist_dip <- dist((matrix_dip))
euc_dist_tetra <- dist((matrix_tetra))

euc_dist_dip_wardd2 <- (hclust(d= euc_dist_dip, method = "ward.D2"))
euc_dist_tetra_wardd2 <- (hclust(d= euc_dist_tetra, method = "ward.D2"))

# create tanglegram
dend1 <- as.dendrogram(euc_dist_dip_wardd2) 
dend2 <- as.dendrogram(euc_dist_tetra_wardd2)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = "Diploid vs. tetraploid, WardD2 using euclidian distance",
           cex_main = 1.5,
           sub = paste("entanglement =", round(entanglement(dend_list), 2)),
           cex_sub = 1,
           lab.cex = 0.5
)

# tanglegram comparing nei and euclidian distance
dend1 <- as.dendrogram(euc_dist_tetra_wardd2) 
dend2 <- as.dendrogram(nei_dist_tetra_wardd2)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = "WardD2 using euclidian vs. Nei's distance",
           cex_main = 1.5,
           sub = paste("entanglement =", round(entanglement(dend_list), 2)),
           cex_sub = 1,
           lab.cex = 0.5
)


# hist(nei_dist_tetra)
# hist(euc_dist_tetra)

# average clustering
euc_dist_average <- (hclust(d= euc_dist_tetra, method = "average"))
nei_dist_average <- (hclust(d= nei_dist_tetra, method = "average"))

dend1 <- as.dendrogram(nei_dist_tetra_wardd2) 
dend2 <- as.dendrogram(nei_dist_average)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = "Average vs. Ward.D2 clustering using Nei's distance",
           cex_main = 1.5,
           sub = paste("entanglement =", round(entanglement(dend_list), 2)),
           cex_sub = 1,
           lab.cex = 0.5
)


dend1 <- as.dendrogram(euc_dist_tetra_wardd2) 
dend2 <- as.dendrogram(euc_dist_average)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = "Average vs. Ward.D2 clustering using euclidian distance",
           cex_main = 1.5,
           sub = paste("entanglement =", round(entanglement(dend_list), 2)),
           cex_sub = 1,
           lab.cex = 0.5
)

# euclidian distance agrees with ADMIXTURE to a higher degree when looking at the subset > 0.8
dend <- as.dendrogram(euc_dist_tetra_wardd2) 
#dend <- as.dendrogram(nei_dist_tetra_wardd2) 
labels_colors(dend) <- admixture$subpopulation[order.dendrogram(dend)]

labels_cex(dend) <- 0.5
# for data not subset
png(filename = "dendrogram_with_subpopulations.png", width = 4000, height = 1000, units = "px", res = 200)
dend %>% set("branches_k_color", value = c("red", "blue", "green", "black", "cyan"), k = 5) %>% 
  #set("labels_col", value = c("skyblue", "orange", "grey", "red", "green"), k = 5) %>% # colors clusters, not according to subpopulations. Gonna have to use default colors
  plot
dev.off()

# for subset data
png(filename = "dendrogram_with_subpopulations>0.8.png", width = 2000, height = 1000, units = "px", res = 200)
dend %>% set("branches_k_color", value = c("blue", "red", "green", "cyan", "black"), k = 5) %>% plot
dev.off()




# PCA with subpopulations and admixed individuals colored. 
# Scaling isn't necessary since all the variables are the same type of data. With scaling, some individuals cluster extremely far from the rest
pca = prcomp((matrix_tetra), scale. = F)

pc1 <- autoplot(pca, label = T,label.repel=T, data = admixture, colour = 'subpopulation')+theme()
pc2 <- autoplot(pca, x=2, y=3, label = T,label.repel=T, data = admixture, colour = 'subpopulation')
#pc3 <- autoplot(pca, x=3, y=4, label = T,label.repel=T, data = admixture, colour = 'subpopulation')
#pc4 <- autoplot(pca, x=4, y=5, label = T,label.repel=T, data = admixture, colour = 'subpopulation')

png(filename = "PCA_with_subpopulations.png", width = 3000, height = 1500, units = "px", res = 200)
ggarrange(pc1,pc2)
dev.off()

# heatmap
# better use a subset of the breeds, 16GB of RAM isn't enough otherwise
library("pheatmap") 
png(filename = "heatmap.png", width = 1200, height = 800, units = "px")
pheatmap(matrix_tetra, color = rev(hcl.colors(50, "Sunset")), cex = 0.8)
dev.off()

# distance matrix
png(filename = "distance_matrix.png", width = 2000, height = 2000, units = "px", res = 150)
fviz_dist(dist(matrix_tetra), lab_size = 5)
dev.off()