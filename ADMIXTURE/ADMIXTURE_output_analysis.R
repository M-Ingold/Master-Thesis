library(gridExtra)

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

png(filename = "ADMIXTURE_3_clusters_all_SNPs.png", width = 2000, height = 800, units = "px", res = 100)
par(mar = c(7, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted1)), col=rainbow(3),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 1.5, cex.axis = 1.5, las=2,
        cex = 0.6, cex.main = 2,
        main = "ADMIXTURE results of 3 clusters without MAF filtering"
        )
dev.off()

tbl3=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.5.Q")
rownames(tbl3) <- Samples$VARIETY
tbl_sorted3 <- tbl3[order(tbl3$V1,tbl3$V4,tbl3$V2,tbl3$V3,tbl3$V5),]

png(filename = "ADMIXTURE_5_clusters_all_SNPs.png", width = 2000, height = 800, units = "px", res = 100)
par(mar = c(7, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted3)), col=rainbow(5),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 1.5, cex.axis = 1.5, las=2,
        cex = 0.6, cex.main = 2,
        main = "ADMIXTURE results of 5 clusters without MAF filtering"
)
dev.off()


tbl2=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.5.Q")
rownames(tbl2) <- Samples$VARIETY
tbl_sorted2 <- tbl2[order(tbl2$V1,tbl2$V4,tbl2$V3,tbl2$V5,tbl2$V4),]

png(filename = "ADMIXTURE_5_clusters_MAF.png", width = 2000, height = 800, units = "px", res = 100)
par(mar = c(7, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted2)), col=rainbow(5),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 1.5, cex.axis = 1.5, las=2,
        cex = 0.6, cex.main = 2,
        main = "ADMIXTURE results of 5 clusters with MAF > 0.05"

)
dev.off()

tbl4=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.3.Q")
rownames(tbl4) <- Samples$VARIETY
tbl_sorted4 <- tbl4[order(tbl4$V1,tbl4$V3,tbl4$V2),]

png(filename = "ADMIXTURE_3_clusters_MAF.png", width = 2000, height = 800, units = "px", res = 100)
par(mar = c(7, 8, 4, 2))
barplot(t(as.matrix(tbl_sorted4)), col=rainbow(3),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 1.5, cex.axis = 1.5, las=2,
        cex = 0.6, cex.main = 2,   
        main = "ADMIXTURE results of 3 clusters with MAF > 0.05"
        
)
dev.off()


test=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.4.Q")
rownames(test) <- Samples$VARIETY
test_sorted <- test[order(test$V1,test$V4,test$V3),]

par(mar = c(7, 8, 4, 2))
barplot(t(as.matrix(test_sorted)), col=rainbow(4),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 1.5, cex.axis = 1.5, las=2,
        cex = 0.6, cex.main = 2,   
        main = "ADMIXTURE results of 4 clusters with MAF > 0.05"
        
)


# pophelper
library(pophelper)

sfiles <- list.files(path = "../scripts/ADMIXTURE/", pattern = "freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.*.Q", full.names=T)
slist <- sortQ(readQ(files=sfiles))
tr1 <- tabulateQ(qlist=slist)

png(filename = "ADMIXTURE_1_to_10_clusters_all_SNPs.png", width = 2000, height = 2000, units = "px")
p1 <- plotQ(slist,imgoutput="join",returnplot=T,exportplot=F,basesize=20, sortind = "all", sharedindlab = F)
grid.arrange(p1$plot[[1]])
dev.off()



# plot maximum subpopulation membership in descending order
tbl2$max <- pmax(tbl2$V1,tbl2$V2,tbl2$V3,tbl2$V4,tbl2$V5)
tbl2 <- tbl2[order(tbl2$max, decreasing = T),]$max
plot(tbl2, ylab = "Maximum subpopulation membership", xlab = '')
abline(h=0.8, lwd = 2, col = 'red')
