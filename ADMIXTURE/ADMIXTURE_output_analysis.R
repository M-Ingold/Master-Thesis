

#copy - pasted values for CV error depending on number of clusters
K <- c(0.19391,0.18748,0.18594,0.18442,0.18390,0.18397,0.18267,0.18567,0.18342,0.18500)

png(filename = "ADMIXTURE_cluster_number_all_SNPs.png")
plot(K,xlab="Number of Clusters", ylab="Cross-Validation error", main="ADMIXTURE CV error by number of clusters")
dev.off()

#copy - pasted values for CV error depending on number of clusters
K2 <- c(0.37532,0.36856,0.36282,0.36353,0.36019,0.35873,0.35604,0.35374,0.35267,0.35269)

png(filename = "ADMIXTURE_cluster_number_MAF.png")
plot(K2,xlab="Number of Clusters", ylab="Cross-Validation error", main="ADMIXTURE CV error by number of clusters")
dev.off()



tbl1=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.3.Q")
rownames(tbl1) <- Samples$VARIETY
tbl_sorted1 <- tbl1[order(tbl1$V1,tbl1$V3,tbl1$V2),]

png(filename = "ADMIXTURE_3_clusters_all_SNPs.png", width = 2000, height = 800, units = "px")
par(mar = c(5, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted1)), col=rainbow(3),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 2, cex.axis = 1.5, las=2,
        cex = 0.5,
        main = "ADMIXTURE results of 3 clusters without MAF filtering"
        )
dev.off()

tbl3=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.5.Q")
rownames(tbl3) <- Samples$VARIETY
tbl_sorted3 <- tbl3[order(tbl3$V1,tbl3$V4,tbl3$V2,tbl3$V3,tbl3$V5),]

png(filename = "ADMIXTURE_5_clusters_all_SNPs.png", width = 2000, height = 800, units = "px")
par(mar = c(5, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted3)), col=rainbow(5),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 2, cex.axis = 1.5, las=2,
        cex = 0.5,
        main = "ADMIXTURE results of 5 clusters without MAF filtering"
)
dev.off()


tbl2=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.5.Q")
rownames(tbl2) <- Samples$VARIETY
tbl_sorted2 <- tbl2[order(tbl2$V1,tbl2$V4,tbl2$V2,tbl2$V3,tbl2$V5),]

png(filename = "ADMIXTURE_5_clusters_MAF.png", width = 2000, height = 800, units = "px")
par(mar = c(5, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted2)), col=rainbow(5),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 2, cex.axis = 1.5, las=2,
        cex = 0.5,       
        main = "ADMIXTURE results of 5 clusters with MAF > 0.05"

)
dev.off()

tbl4=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.3.Q")
rownames(tbl4) <- Samples$VARIETY
tbl_sorted4 <- tbl4[order(tbl4$V1,tbl4$V3,tbl4$V2),]

png(filename = "ADMIXTURE_3_clusters_MAF.png", width = 2000, height = 800, units = "px")
par(mar = c(5, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted4)), col=rainbow(3),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 2, cex.axis = 1.5, las=2,
        cex = 0.5,       
        main = "ADMIXTURE results of 3 clusters with MAF > 0.05"
        
)
dev.off()

