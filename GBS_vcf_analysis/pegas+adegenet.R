library(pegas)
#library(snpStats)
library(adegenet)
library(ggplot2)
#library(GGtools)
library(poppr)

#x <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf', to = 50000)
x <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf', to = 50000)
#x <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf', to = 50000)
#x <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf', to = 50000)

rownames(x) <- Samples$VARIETY

# assign population. Not useful this way for AMOVA
newx <- merge(x, admixture, by = 0, all = F)
newx$V1 <- NULL
newx$V2 <- NULL
newx$V3 <- NULL
newx$V4 <- NULL
newx$V5 <- NULL
rownames(newx) <- sort(rownames(admixture))
newx$Row.names <- NULL
newx$subpopulation

X <- as.loci(newx, col.pop = length(newx), ploidy = 4
             )
X <- as.loci(x, ploidy = 4
            )
#convert to adegenet genind
y <- loci2genind(X, ploidy = 4)
#y <- loci2genind(X, ploidy = 2)

grp <- find.clusters(y, scale = F)
xvalDapc(y, grp$grp)
dapc <- dapc(y, grp$grp, scale = F, center = T, n.pca = 200, n.da = 2)
scatter(dapc)
compoplot(dapc, posi="bottomright",
          txt.leg=paste("Cluster", 1:5), lab="",
          ncol=2, xlab="individuals", col=funky(5))

# subpopulations need to be added this way as strata to be able to perform AMOVA
strata(y) <- admixture[order(row.names(admixture)),]$subpopulation
#strata(y) <- admixture[order(row.names(admixture)),]

#convert to genpop
#z <- genind2genpop(y)

#convert to genclone
gc <- poppr::as.genclone(y)

amova.result <- poppr.amova(gc, ~subpopulation)
amova.result
amova.test <- randtest(amova.result)
plot(amova.test)
amova.test
# to-do: reshuffle as test

# probably bad way of doing this
# hamming distance takes very long
#dist <- dist.hamming(X)
# euclidean <- dist(X)
# res.amova <- amova(euclidean ~ population, data = X)
# res.amova

# # test
# Y <- sample(c("A", "G"), size = 45, replace = TRUE)
# Y <- matrix(Y, 9, 5)
# Y
# Yloc <- as.loci(as.data.frame(Y))
# Yloc$population <- gl(3, 3, labels = paste0("pop", 1:3))
# Yg <- loci2genind(Yloc, ploidy = 1)
# dy <- dist.hamming(Yloc)
# dy
# res.amova <- amova(dy ~ population, data = Yloc)
# 
# data(Aeut)
# strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
# agc <- as.genclone(Aeut)
# agc
# amova.result <- poppr.amova(agc, ~Pop/Subpop)
# amova.result
# amova.test <- randtest(amova.result) # Test for significance
# plot(amova.test)

#subset populations
pops <- seppop(y)

# HWE of the populations
png(filename = "HWE_of_populations_tetraploid.png", width=1000, height = 1500, res = 120)
n = 0
par(mfrow=c(3,2))
for (pop in pops) {
  n <- n+1
  hwp <- hw.test(pop, B=0)
  hist(hwp[,3], breaks = 100, main = paste("chi²-test for HWE of population ", n), xlab = "p-value")
}
dev.off()

# HWE of the populations
png(filename = "MAFs_of_populations_tetraploid.png", width=1000, height = 1500, res = 120)
n = 0
par(mfrow=c(3,2))
for (pop in pops) {
  n <- n+1
  maf <- minorAllele(pop)
  hist(maf, breaks = 100, main = paste("MAFs of population", n), xlab = "MAF", xlim = c(0,0.5))
}
dev.off()

hist(minorAllele(y), breaks = 100, xlim = c(0,0.5))

#expected vs observed heterozygosity of the populations
png(filename = "heterozygosity_of_populations_diploid.png", width=1500, height = 2000)
n = 0
par(mfrow=c(5,2))
for (pop in pops) {
  n <- n+1
  sum <- summary(pop)
  hist(sum$Hobs, main = paste("Histogram of observed heterozygosity of population", n), cex = 30)
  abline(v = mean(sum$Hobs), col="red", lwd=3, lty=2)
  mtext(round(mean(sum$Hobs), 2),side = 3)

  hist(sum$Hexp, main = paste("Histogram of expected heterozygosityof population", n), cex = 30)
  abline(v = mean(sum$Hexp), col="red", lwd=3, lty=2)
  mtext(round(mean(sum$Hexp), 2),side = 3)
  }
dev.off()


# #expected vs observed heterozygosity might be interesting
# sumy <- summary(y)
#  
# png(filename = "heterozygosity.png", width=1000, height = 500)
# par(mfrow=c(1,2))
# hist(sumy$Hobs, main = "Histogram of observed heterozygosity")
# abline(v = mean(sumy$Hobs), col="red", lwd=3, lty=2)
# mtext(round(mean(sumy$Hobs), 2),side = 3)
# hist(sumy$Hexp, main = "Histogram of expected heterozygosity")
# abline(v = mean(sumy$Hexp), col="red", lwd=3, lty=2)
# mtext(round(mean(sumy$Hexp), 2),side = 3)
# 
# dev.off()
# 
hwp <- hw.test(y, B=0)

png(filename = "chi²-pval-hist.png", width=1000, height = 1000)
hist(hwp[,3], breaks = 100, main = "chi²-test for HWE", xlab = "p-value")
dev.off()

# only diploid data
# Fst(x)
# Rst(X)


#snpStats: hard to import using GGtools
#ss <- vcf2sm('../data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf')



# plot HWE onto genome
# probably no interesting pattern
# hw_df <- as.data.frame(hwp)
# colnames(hw_df) <- c("chisq","df", "Pr_chisq")
# 
# vcf <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf")
# hw_df$chr <- vcf@fix[,1]
# hw_df$pos <- vcf@fix[,2]
# 
# p <- ggplot(hw_df, aes(pos,Pr_chisq)) + geom_point()
# 
# png(filename = "hwe_plotted.png", width=1000, height = 1000)
# p + facet_wrap(vars(chr))
# dev.off()
# 
# # some plots seem to be inverted compared to ggplot?
# hwe1 <- subset(hw_df, chr=="chr01")
# plot(hwe1$pos, hwe1$Pr_chisq)
# 
# hwe2 <- subset(hw_df, chr=="chr02")
# plot(hwe2$pos, hwe2$Pr_chisq)
# 
# hwe11 <- subset(hw_df, chr=="chr11")
# plot(hwe11$pos, hwe11$Pr_chisq)
# 
# 

# population strata needed?
# amova <- poppr.amova(y)

# only diploids
#fs <- filter_stats(y, plot = T)

# # Index of Association
# ia(y)
# 
# #long runtime
# res <- pair.ia(y)
# plot(res, low = "black", high = "green", index = "Ia")

# needs snplight
# winia <- win.ia(y)
# plot(winia, low = "black", high = "green", index = "Ia")


# loctable <- locus_table(y)
# #Hexp
# hist(loctable[,3])
# 
# #Evenness?
# hist(loctable[,4])
# 
# #1-D?
# hist(loctable[,2])

# summary statistics, sample for p-value takes a while
summary <- poppr(y
                 #,sample = 999 # for p-values
                 )

# Probability a genotype is derived from sexual reproduction
# only for haploid or diploid data
#psex(y)

# Fst by hand
Ht <- summary$Hexp[6]
Hsub <- summary[-6,c(2,10)] # -6 or -7 depending on if admixed individuals are counted as a subpopulation
Hs <- sum(Hsub[1]*Hsub[2])/sum(Hsub[1])

Fst <- (Ht - Hs)/Ht
# 4.7% of variation is explained by the subpopulations when MAF > 0.05
# 5% of variation is explained by the subpopulations when MAF wasn't filtered
# 3.2% of variation is explained by the subpopulations when MAF MAF > 0.05 for diploidized data
# 3.2% of variation is explained by the subpopulations when MAF wasn't filtered for diploidized data

#to-do: Fis

# inbreeding coefficient
# everything is loaded fresh to not exclude any samples as before

#x <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf', to = 50000)
#x <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf', to = 50000)
x <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf', to = 50000)
#x <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf', to = 50000)


X <- as.loci(x)

y <- loci2genind(X, ploidy = 4)

# inb <- inbreeding(y, res.type = "estimate")
# 
# Fmean=sapply(inb, mean)
# hist(Fmean, col="orange", xlab="mean value of F",
# main="Distribution of mean F across individuals")

# calculate 
sum <- summary(y)


HT <- mean(sum$Hexp) #without subpopulations equivalent to Hs
HI <- mean(sum$Hobs)
FIT <- (HT - HI)/HT

FIT2 <- 1-(HI/HT)



library(hierfstat)

# I don't get why this doesn't work
# Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  : arguments imply differing number of rows: 75, 74
basic.stats(y)
hfstat <- genind2hierfstat(y)
pairwise.neifst(hfstat,diploid=TRUE)

# #test
# data(nancycats)
# genind2hierfstat(nancycats)
# basic.stats(nancycats)

#directly read VCF, populations need to be added later
hfstat <- read.VCF('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf', convert.chr = F)
# how do i need to process the bed.matrix to use this method?
basic.stats(hfstat@snps)

library(SNPRelate)

snpgdsVCF2GDS('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf', "test.gds")
genofile <- openfn.gds("test.gds")
read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "genotype"))

#seqSetFilter(genofile, sample.id=samp.id[c(2,4,6)])

pca <- snpgdsPCA(genofile, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
