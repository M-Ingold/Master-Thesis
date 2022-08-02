library(gridExtra)

# ADMIXTURE results of MAF-filtered and unfiltered VCFs are visualized

# copy - pasted values for CV error depending on number of clusters
K <- c(0.19391,0.18748,0.18594,0.18442,0.18390,0.18397,0.18267,0.18567,0.18342,0.18500)
K_maf <- c(0.37532,0.36856,0.36282,0.36353,0.36019,0.35873,0.35604,0.35374,0.35267,0.35269)

png(filename = "ADMIXTURE_CV_error.png", width = 1000, height = 500, res=100)
par(mfrow = c(1,2))
plot(K_maf,xlab="Number of Clusters", ylab="Cross-Validation error", main="With MAF > 0.5 filter")
plot(K,xlab="Number of Clusters", ylab="Cross-Validation error", main="Without MAF filter")
dev.off()


# 4 barplots of ADMIXTURE results are saved in one plot. Optional code to save plots separately is commented out

png(filename = "ADMIXTURE_4_plots.png", width = 2000, height = 2400, units = "px", res = 100)
par(mfrow = c(4,1))
par(mar = c(0, 8, 4, 0))

#png(filename = "ADMIXTURE_3_clusters_all_SNPs.png", width = 2000, height = 800, units = "px", res = 100)

tbl4=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.3.Q")
#rownames(tbl4) <- Samples$VARIETY
tbl_sorted4 <- tbl4[order(tbl4$V1,tbl4$V3,tbl4$V2),]

#png(filename = "ADMIXTURE_3_clusters_MAF.png", width = 2000, height = 800, units = "px", res = 100)
barplot(t(as.matrix(tbl_sorted4)), col=rainbow(3),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 3, cex.axis = 1.5,
        cex = 1, cex.main = 3,
        main = "ADMIXTURE results of 3 clusters with MAF > 0.05",
        xaxt='n'
        
)


tbl2=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.5.Q")
#rownames(tbl2) <- Samples$VARIETY
tbl_sorted2 <- tbl2[order(tbl2$V1,tbl2$V4,tbl2$V2,tbl2$V3,tbl2$V4),]

#png(filename = "ADMIXTURE_5_clusters_MAF.png", width = 2000, height = 800, units = "px", res = 100)
par(mar = c(0, 8, 4, 0))
barplot(t(as.matrix(tbl_sorted2)), col=rainbow(5),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 3, cex.axis = 1.5,
        cex = 1, cex.main = 3,
        main = "ADMIXTURE results of 5 clusters with MAF > 0.05",
        xaxt='n'

)
#dev.off()

tbl1=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.3.Q")
#rownames(tbl1) <- Samples$VARIETY
tbl_sorted1 <- tbl1[order(tbl1$V1,tbl1$V3,tbl1$V2),]

barplot(t(as.matrix(tbl_sorted1)), col=rainbow(3),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 3, cex.axis = 1.5,
        cex = 1, cex.main = 3,
        main = "ADMIXTURE results of 3 clusters without MAF filtering",
        xaxt='n'
        )
#dev.off()

tbl3=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.5.Q")
rownames(tbl3) <- Samples$VARIETY
tbl_sorted3 <- tbl3[order(tbl3$V1,tbl3$V4,tbl3$V2,tbl3$V3,tbl3$V5),]

#png(filename = "ADMIXTURE_5_clusters_all_SNPs.png", width = 2000, height = 800, units = "px", res = 100)
barplot(t(as.matrix(tbl_sorted3)), col=rainbow(5),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 3, cex.axis = 1.5,
        cex = 1, cex.main = 3,
        main = "ADMIXTURE results of 5 clusters without MAF filtering",
        xaxt='n'
)
#dev.off()



dev.off()

# MAF-filtered VCF results for four populations is used as separate plot

test=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf.4.Q")
#rownames(test) <- Samples$VARIETY
test_sorted <- test[order(test$V1,test$V4,test$V3),]

png(filename = "ADMIXTURE_4_clusters_MAF.png", width = 2000, height = 800, units = "px", res = 100)
par(mar = c(1, 8, 4, 0))
barplot(t(as.matrix(test_sorted)), col=rainbow(4),
        #xlab="Individual #", 
        ylab="Proportional Membership (Q)", border=NA,
        cex.lab = 3, cex.axis = 1.5,
        cex = 1, cex.main = 3,
        las=2,
        main = "ADMIXTURE results of 4 clusters with MAF > 0.05"
        
)
dev.off()


# pophelper is a nice package to visualize results for all numbers of K, not used in the master thesis

library(pophelper)

sfiles <- list.files(path = "../scripts/ADMIXTURE/", pattern = "freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.*.Q", full.names=T)
slist <- sortQ(readQ(files=sfiles))
tr1 <- tabulateQ(qlist=slist)

png(filename = "ADMIXTURE_1_to_10_clusters_all_SNPs.png", width = 2000, height = 2000, units = "px")
p1 <- plotQ(slist,imgoutput="join",returnplot=T,exportplot=F,basesize=20, sortind = "all", sharedindlab = F)
grid.arrange(p1$plot[[1]])
dev.off()



# plot maximum subpopulation membership in descending order

tbl3$max <- pmax(tbl3$V1,tbl3$V2,tbl3$V3,tbl3$V4,tbl3$V5)
rownames(tbl3) <- Samples$VARIETY
max_membership <- tbl3[order(tbl3$max, decreasing = T),]$max
png(filename = "max_subpopulation_membership", width = 500, height = 500, res = 100)
plot(max_membership, pch=19, ylab = "Maximum subpopulation membership", xlab = '')
abline(h=0.8, lwd = 2, col = 'red')
dev.off()