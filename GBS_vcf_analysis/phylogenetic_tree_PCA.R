library(ape)
library(vcfR)

vcfTetra <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf")
dosagetetra <- extract.gt(vcfTetra, as.numeric = F)

dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/0", 0)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/1", 1)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/1/1", 2)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/1/1/1", 3)
dosagetetra <- replace(dosagetetra, dosagetetra == "1/1/1/1", 4)
dosagetetra <- as.numeric(na.omit(dosagetetra))

testmatrix <- matrix(unlist(dosagetetra), nrow = 5187, ncol = 261)

pca = prcomp(t(testmatrix), scale. = T)
vars = 100*pca$sdev / sum(pca$sdev)
plot(pca$x[,1],pca$x[,2],xlab = sprintf("PCA1 (%.1f%%)", pca.var.per[1]), ylab = sprintf("PCA2 (%.1f%%)", pca.var.per[2]))
plot(pca$x[,2],pca$x[,3],xlab = sprintf("PCA2 (%.1f%%)", pca.var.per[2]), ylab = sprintf("PCA3 (%.1f%%)", pca.var.per[3]))

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")



# Genomic distance
# scaling not necessary because no different units are measured?
ss <- sample(1:261, 20)
dist <- dist(t(testmatrix[,ss]))
dist <- dist(scale(t(testmatrix)))
dist <- dist(t(testmatrix))

# phylogram
tre <- nj(dist)
par(xpd=TRUE)
plot(tre, type="phylogram")
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))

# distance matrix
fviz_dist(dist)

# Number of clusters zero, or 2 or 3 without transposition?
fviz_nbclust(t(testmatrix), kmeans, method = "wss", k.max = 100)

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
pam.res <- pam(t(testmatrix), 2, medoids = c(155, 33)) 
print(pam.res)

# clusters not in line with kmeans
# due to clustering along other variables than PC1?
fviz_cluster(pam.res, palette = c("#00AFBB", "#FC4E07"), # color palette 
             ellipse.type = "t", # Concentration ellipse 
             repel = TRUE, # Avoid label overplotting (slow) 
             ggtheme = theme_classic() 
)

