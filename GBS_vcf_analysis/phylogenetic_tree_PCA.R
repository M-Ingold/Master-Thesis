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

# to-do: test without admixed individuals
#        save plots

vcfTetra <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf")
vcfDiploid <- read.vcfR("../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf")

dosagetetra <- extract.gt(vcfTetra, as.numeric = F)

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


# get diploid dosages. Alternatively one could diploidize the tetraploid VCF
dosagedip <- extract.gt(vcfDiploid, as.numeric = F)

dosagedip <- replace(dosagedip, dosagedip == "0/0", 0)
dosagedip <- replace(dosagedip, dosagedip == "0/1", 1)
dosagedip <- replace(dosagedip, dosagedip == "1/1", 2)
dosagedip[is.na(dosagedip)] <- 0
dosagedip <- as.numeric(dosagedip)

matrix_dip <- t(matrix(unlist(dosagedip), nrow = length(dosagedip)/261, ncol = 261))


# import potato genotype data
GBS02 <- read_excel("../data/GBS2/GBS02_sample-IDs_14.04.21.xlsx", skip = 6)
GBS02 <- rename(GBS02, Identifier = "Sample identifier")
GBS01 <- read_excel("../data/GBS1/22-01-2021_Zuordnung_genotypes_SB.aktualisiert.xlsx")

# poorly merge datasets, maybe do it better later. Too bad!
intersect <- intersect(names(GBS01), names(GBS02))
GBS <- merge(GBS01, GBS02, by=intersect, all = T)

df <- data.frame(Identifier=gsub("Sample_", "",colnames(vcfTetra@gt)[-1]))
Samples <- merge(df, GBS, by="Identifier")

# add sample code to matrix to later substitute proper breed name
rownames(matrix_tetra) <- Samples$VARIETY
rownames(matrix_dip) <- Samples$VARIETY

# tanglegrams of diploid VCF compared to tetraploid VCF
library(poppr)
res.dist_dip <- nei.dist((matrix_dip))
res.dist_tetra <- nei.dist((matrix_tetra))

res.hc_dip_wardd2 <- hclust(d= res.dist_dip, method = "ward.D2")
res.hc_tetra_wardd2 <- hclust(d= res.dist_tetra, method = "ward.D2")

# Create two dendrograms 
dend1 <- as.dendrogram(res.hc_dip_wardd2) 
dend2 <- as.dendrogram(res.hc_tetra_wardd2)

# Create a list to hold dendrograms 
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2)) 
)

fviz_dend(res.hc_dip_wardd2, k = 5, # Cut in groups 
          cex = 0.5, # label size 
          #k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
          color_labels_by_k = TRUE, # color labels by groups 
          rect = F, # Add rectangle around groups 
          #type = "circular"
)

res.hc_dip_average <- hclust(d= res.dist_dip, method = "average")
res.hc_tetra_average <- hclust(d= res.dist_tetra, method = "average")

dend3 <- as.dendrogram(res.hc_dip_average) 
dend4 <- as.dendrogram(res.hc_tetra_average)

dend_list <- dendlist(dend3, dend4)

tanglegram(dend3, dend4, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2)) 
)




# subset Samples by ADMIXTURE population membership > 90%
admixture=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.5.Q")
rownames(admixture) <- Samples$VARIETY
admixture <- subset(admixture, V1 > 0.9|V2 > 0.9|V3 > 0.9|V4 > 0.9|V5 > 0.9)

# create dataframe row containing population number
admixture$subpopulation <- NA
for (breed in 1:nrow(admixture)){
  if (admixture[breed,]$V1 > 0.9){
    admixture[breed,]$subpopulation <- "1"
  }
  if (admixture[breed,]$V2 > 0.9){
    admixture[breed,]$subpopulation <- "2"
  }
  if (admixture[breed,]$V3 > 0.9){
    admixture[breed,]$subpopulation <- "3"
  }
  if (admixture[breed,]$V4 > 0.9){
    admixture[breed,]$subpopulation <- "4"
  }
  if (admixture[breed,]$V5 > 0.9){
    admixture[breed,]$subpopulation <- "5"
  }
}

#admixture_subpops <- admixture$subpopulation
#rownames(admixture_subpops) <- rownames(admixture)
#write.csv(admixture, file = "75_non-admixed_breeds.csv")

# subset for non-admixed breeds
matrix_tetra <- matrix_tetra[rownames(admixture), ]
matrix_dip <- matrix_dip[rownames(admixture), ]

# # subset for non-admixed breeds with known maturity etc
# matrix_tetra <- matrix_tetra[x$Row.names, ]
# matrix_dip <- matrix_dip[rownames(x), ]

# after subsetting, some rows contain all the same number. These have to be removed, in this case using the variance.
matrix_tetra <- (matrix_tetra)[ , which(apply((matrix_tetra), 2, var) != 0)]
matrix_dip <- (matrix_dip)[ , which(apply((matrix_dip), 2, var) != 0)]


pca = prcomp((matrix_tetra), scale. = T)
pca_dip = prcomp((matrix_dip), scale. = T)

vars = 100*pca$sdev / sum(pca$sdev)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
plot(pca$x[,1],pca$x[,2],xlab = sprintf("PC1 (%.1f%%)", pca.var.per[1]), ylab = sprintf("PC2 (%.1f%%)", pca.var.per[2]))
plot(pca$x[,2],pca$x[,3],xlab = sprintf("PC2 (%.1f%%)", pca.var.per[2]), ylab = sprintf("PC3 (%.1f%%)", pca.var.per[3]))
plot(pca$x[,3],pca$x[,4],xlab = sprintf("PC3 (%.1f%%)", pca.var.per[3]), ylab = sprintf("PC4 (%.1f%%)", pca.var.per[4]))

png(filename = "pca_scree_plot.png", width = 500, height = 500, units = "px")
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
dev.off()

pc1 <- autoplot(pca, label = T,label.repel=T, data = admixture, colour = 'subpopulation')
pc2 <- autoplot(pca, x=2, y=3, label = T,label.repel=T, data = admixture, colour = 'subpopulation')
#autoplot(pca, x=3, y=4, label = T,label.repel=T, data = admixture, colour = 'subpopulation')

png(filename = "PCA.png", width = 1000, height = 500, units = "px")
ggarrange(pc1,pc2)
dev.off()

# pc1 <- autoplot(pca, label = T,label.repel=T, data = x, colour = 'COUNTRY')
# pc2 <- autoplot(pca, x=2, y=3, label = T,label.repel=T, data = x, colour = 'COUNTRY')
# ggarrange(pc1,pc2)
# 
# pc1 <- autoplot(pca, label = T,label.repel=T, data = x, colour = 'heat_sensitivity')
# pc2 <- autoplot(pca, x=2, y=3, label = T,label.repel=T, data = x, colour = 'heat_sensitivity')
# ggarrange(pc1,pc2)
# 
# pc1_dip <- autoplot(pca_dip, label = T,label.repel=T, data = admixture, colour = 'subpopulation')
# pc2_dip <- autoplot(pca_dip, x=2, y=3, label = T,label.repel=T, data = admixture, colour = 'subpopulation')
# ggarrange(pc1, pc2, pc1_dip,pc2_dip)

# tanglegram of non-admixed breeds
res.dist_dip <- nei.dist((matrix_dip))
res.dist_tetra <- nei.dist((matrix_tetra))

res.hc_dip_wardd2 <- hclust(d= res.dist_dip, method = "ward.D2")
res.hc_tetra_wardd2 <- hclust(d= res.dist_tetra, method = "ward.D2")

# Create two dendrograms 
dend1 <- as.dendrogram(res.hc_dip_wardd2) 
dend2 <- as.dendrogram(res.hc_tetra_wardd2)

# Create a list to hold dendrograms 
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2)) 
)

fviz_dend(res.hc_dip_wardd2, k = 5, # Cut in groups 
          cex = 0.7, # label size 
          #k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
          color_labels_by_k = TRUE, # color labels by groups 
          rect = F, # Add rectangle around groups 
          #type = "circular"
)

res.hc_dip_average <- hclust(d= res.dist_dip, method = "average")
res.hc_tetra_average <- hclust(d= res.dist_tetra, method = "average")

dend3 <- as.dendrogram(res.hc_dip_average) 
dend4 <- as.dendrogram(res.hc_tetra_average)

dend_list <- dendlist(dend3, dend4)

png(filename = "tanglegram_average_diploid_tetraploid.png", width = 1200, height = 800, units = "px")
tanglegram(dend3, dend4, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2)) 
)
dev.off()


# Genomic distance
# scaling not necessary because no different units are measured?
#ss <- sample(1:261, 20)
#dist <- dist(t(matrix_tetra[,ss]))
#dist <- dist(scale(t(matrix_tetra)))
dist <- nei.dist((matrix_tetra))
res.dist <- nei.dist((matrix_tetra))

# phylogram
tre <- nj(dist)
par(xpd=TRUE)
plot(tre, type="phylogram")
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))

# distance matrix
png(filename = "distance_matrix.png", width = 1000, height = 1000, units = "px")
fviz_dist(dist)
dev.off()

# Number of clusters zero, or 2 or 3 without transposition? 5 without admixed breeds, maybe 4?
png(filename = "k-means_cluster_number.png", width = 700, height = 500, units = "px")
fviz_nbclust((matrix_tetra), kmeans, method = "wss", k.max = 10)
dev.off()

km.res <- kmeans((matrix_tetra), 2, nstart = 25)
km.res4 <- kmeans((matrix_tetra), 4, nstart = 25)
km.res5 <- kmeans((matrix_tetra), 5, nstart = 25)

png(filename = "k-means_cluster_plot.png", width = 800, height = 500, units = "px")
fviz_cluster(km.res4, data = (matrix_tetra), #palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
             ellipse.type = "euclid", # Concentration ellipse 
             star.plot = TRUE, # Add segments from centroids to items 
             repel = TRUE, # Avoid label overplotting (slow) 
             main = "Kmeans Cluster Plot",
             ggtheme = theme_minimal() 
)
dev.off()

# PAM
# 2 clusters
fviz_nbclust((matrix_tetra), pam, method = "silhouette")
pam.res <- pam((matrix_tetra), 5) 
print(pam.res)

# clusters not in line with kmeans
# due to clustering along other variables than PC1?
fviz_cluster(pam.res, #palette = c("#00AFBB", "#FC4E07"), # color palette 
             ellipse.type = "t", # Concentration ellipse 
             repel = TRUE, # Avoid label overplotting (slow) 
             ggtheme = theme_classic() 
)

# hierarcical clustering
res.hc <- hclust(d= res.dist, method = "ward.D2")

png(filename = "Ward.D2 clustering.png", width = 1000, height = 800, units = "px")
par(mar = c(5, 10, 4, 2))
fviz_dend(res.hc, cex = 0.8, k = 5,
          main = "Ward.D2 clustering")
dev.off()

# Compute cophentic distance 
res.coph <- cophenetic(res.hc)

# Correlation between cophenetic distance and the original distance 
cor(res.dist, res.coph)

res.hc2 <- hclust(res.dist, method = "average") 
cor(res.dist, cophenetic(res.hc2))

png(filename = "Average_clustering.png", width = 1000, height = 800, units = "px")
par(mar = c(5, 10, 4, 2))
fviz_dend(res.hc2, cex = 0.8, k = 5,
          main = "Average clustering")
dev.off()

# Cut tree into groups 
grp <- cutree(res.hc, k = 5) 
head(grp, n= 4)

# Cut in groups and color by groups
fviz_dend(res.hc2, k = 5, # Cut in groups 
          cex = 0.5, # label size 
          #k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
          color_labels_by_k = TRUE, # color labels by groups 
          rect = F, # Add rectangle around groups 
          #type = "circular"
          )

# only works with defined number of groups?
fviz_cluster(list(data = (matrix_tetra), cluster = grp), 
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
             ellipse.type = "convex", # Concentration ellipse 
             repel = TRUE, # Avoid label overplotting (slow) 
             show.clust.cent = FALSE, ggtheme = theme_minimal())

# Agglomerative Nesting (Hierarchical Clustering) 
res.agnes <- agnes(x= (matrix_tetra), # data matrix 
                   stand = FALSE, # Standardize the data 
                   metric = "euclidean", # metric for distance matrix 
                   method = "average" # Linkage method 
                   )

fviz_dend(res.agnes, cex = 0.6, k= 5,type = "circular")

# Create two dendrograms 
dend1 <- as.dendrogram(res.hc) 
dend2 <- as.dendrogram(res.hc2)

# Create a list to hold dendrograms 
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2)

png(filename = "tanglegram_wardD2_average.png", width = 1200, height = 800, units = "px")
tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2)) 
           )
dev.off()

# Cophenetic correlation matrix 
cor.dendlist(dend_list, method = "cophenetic")
# Baker correlation matrix 
cor.dendlist(dend_list, method = "baker")

# unrooted layout
require("igraph") 
fviz_dend(res.hc2, h= 90, k_colors = "npg", type = "phylogenic", repel = TRUE, phylo_layout = "layout.gem")


# Heatmaps
library("RColorBrewer") 
#col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

#heatmap(matrix_tetra, scale = "none", col = col)

library("pheatmap") 
png(filename = "heatmap.png", width = 1200, height = 800, units = "px")
pheatmap(matrix_tetra, color = rev(hcl.colors(50, "Sunset")), cex = 0.8)
dev.off()

#pheatmap(matrix_tetra, color = col, cex = 0.8)

# Cluster validation
library(clustertend)
hopkins((matrix_tetra), n = 260)
# H = 0.378, not clustering great i guess?

fviz_dist(dist((matrix_tetra)), show_labels = T, lab_size = 5)
