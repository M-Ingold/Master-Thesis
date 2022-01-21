library(ape)
library(vcfR)
library(cluster)
library(ggplot2)
library(factoextra)
library(dendextend)
library(readxl)
library(dplyr)

# to-do: test without diploid breeds
#        save plots

vcfTetra <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf")
dosagetetra <- extract.gt(vcfTetra, as.numeric = F)

dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/0", 0)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/1", 1)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/1/1", 2)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/1/1/1", 3)
dosagetetra <- replace(dosagetetra, dosagetetra == "1/1/1/1", 4)
dosagetetra <- as.numeric(na.omit(dosagetetra))

testmatrix <- matrix(unlist(dosagetetra), nrow = 5187, ncol = 261)

# import potato genotype data
GBS02 <- read_excel("/media/rna/CYSTOPHORA/GBS_Analysis/data/GBS2/GBS02_sample-IDs_14.04.21.xlsx", skip = 6)
GBS02 <- rename(GBS02, Identifier = "Sample identifier")
GBS01 <- read_excel("/media/rna/CYSTOPHORA/GBS_Analysis/data/GBS1/22-01-2021_Zuordnung_genotypes_SB.aktualisiert.xlsx")

# poorly merge datasets, maybe do it better later. Good enough for now
intersect <- intersect(names(GBS01), names(GBS02))
GBS <- merge(GBS01, GBS02, by=intersect, all = T)

df <- data.frame(Identifier=gsub("Sample_", "",colnames(vcfTetra@gt)[-1]))
Samples <- merge(df, GBS, by="Identifier")

# add sample code to matrix to later substitute proper breed name
colnames(testmatrix) <- Samples$VARIETY


pca = prcomp(t(testmatrix), scale. = T)
vars = 100*pca$sdev / sum(pca$sdev)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
plot(pca$x[,1],pca$x[,2],xlab = sprintf("PC1 (%.1f%%)", pca.var.per[1]), ylab = sprintf("PC2 (%.1f%%)", pca.var.per[2]))
plot(pca$x[,2],pca$x[,3],xlab = sprintf("PC2 (%.1f%%)", pca.var.per[2]), ylab = sprintf("PC3 (%.1f%%)", pca.var.per[3]))
plot(pca$x[,3],pca$x[,4],xlab = sprintf("PC3 (%.1f%%)", pca.var.per[2]), ylab = sprintf("PC4 (%.1f%%)", pca.var.per[3]))

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")



# Genomic distance
# scaling not necessary because no different units are measured?
ss <- sample(1:261, 20)
dist <- dist(t(testmatrix[,ss]))
dist <- dist(scale(t(testmatrix)))
dist <- dist(t(testmatrix))
res.dist <- dist(t(testmatrix))

# phylogram
tre <- nj(dist)
par(xpd=TRUE)
plot(tre, type="phylogram")
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))

# distance matrix
fviz_dist(dist)

# Number of clusters zero, or 2 or 3 without transposition?
fviz_nbclust(t(testmatrix), kmeans, method = "wss", k.max = 10)

km.res <- kmeans(t(testmatrix), 2, nstart = 25)
km.res3 <- kmeans(t(testmatrix), 3, nstart = 25)

fviz_cluster(km.res3, data = t(testmatrix), palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
             ellipse.type = "euclid", # Concentration ellipse 
             star.plot = TRUE, # Add segments from centroids to items 
             repel = TRUE, # Avoid label overplotting (slow) 
             ggtheme = theme_minimal() 
)

# PAM
# 2 clusters
fviz_nbclust(t(testmatrix), pam, method = "silhouette")
pam.res <- pam(t(testmatrix), 2, #medoids = c(155, 33)
               ) 
print(pam.res)

# clusters not in line with kmeans
# due to clustering along other variables than PC1?
fviz_cluster(pam.res, palette = c("#00AFBB", "#FC4E07"), # color palette 
             ellipse.type = "t", # Concentration ellipse 
             repel = TRUE, # Avoid label overplotting (slow) 
             ggtheme = theme_classic() 
)

# hierarcical clustering
res.hc <- hclust(d= res.dist, method = "ward.D2")
fviz_dend(res.hc, cex = 0.5)

# Compute cophentic distance 
res.coph <- cophenetic(res.hc)

# Correlation between cophenetic distance and the original distance 
cor(res.dist, res.coph)

res.hc2 <- hclust(res.dist, method = "average") 
cor(res.dist, cophenetic(res.hc2))
fviz_dend(res.hc2, cex = 0.5)

# Cut tree into groups 
grp <- cutree(res.hc, h = 90) 
head(grp, n= 4)

# Cut in groups and color by groups
fviz_dend(res.hc2, h = 90, # Cut in groups 
          cex = 0.5, # label size 
          #k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
          color_labels_by_k = TRUE, # color labels by groups 
          rect = F, # Add rectangle around groups 
          #type = "circular"
          )

# only works with defined number of groups?
fviz_cluster(list(data = t(testmatrix), cluster = grp), 
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
             ellipse.type = "convex", # Concentration ellipse 
             repel = TRUE, # Avoid label overplotting (slow) 
             show.clust.cent = FALSE, ggtheme = theme_minimal())

# Agglomerative Nesting (Hierarchical Clustering) 
res.agnes <- agnes(x= t(testmatrix), # data matrix 
                   stand = FALSE, # Standardize the data 
                   metric = "euclidean", # metric for distance matrix 
                   method = "average" # Linkage method 
                   )

fviz_dend(res.agnes, cex = 0.6, h= 90,type = "circular")

# Create two dendrograms 
dend1 <- as.dendrogram (res.hc) 
dend2 <- as.dendrogram (res.hc2)

# Create a list to hold dendrograms 
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2)

tanglegram(dend1, dend2, highlight_distinct_edges = FALSE, # Turn-off dashed lines 
           common_subtrees_color_lines = T, # Turn-off line colors 
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2)) 
           )

# Cophenetic correlation matrix 
cor.dendlist(dend_list, method = "cophenetic")
# Baker correlation matrix 
cor.dendlist(dend_list, method = "baker")

# unrooted layout
require("igraph") 
fviz_dend(res.hc2, h= 90, k_colors = "npg", type = "phylogenic", repel = TRUE, phylo_layout = "layout.gem")


# Heatmaps
library("RColorBrewer") 
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

heatmap(testmatrix, scale = "none", col = col)

library("pheatmap") 
pheatmap(testmatrix, color = rev(hcl.colors(50, "Sunset")), cex = 0.8)
pheatmap(testmatrix, color = col, cex = 0.8)

# Cluster validation
library(clustertend)
hopkins(t(testmatrix), n = 260)
# H = 0.378, not clustering great i guess?

fviz_dist(dist(t(testmatrix)), show_labels = FALSE)
