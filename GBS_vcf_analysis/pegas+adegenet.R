library(pegas)
#library(snpStats)
library(adegenet)
library(ggplot2)
#library(GGtools)
library(poppr)

#x <- read.vcf('../data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf', to = 50000)

x <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf', to = 7000)

X <- as.loci(x)

#convert to adegenet genind
y <- loci2genind(X, ploidy = 4)

#convert to genpop
z <- genind2genpop(y)

#convert to genclone
gc <- poppr::as.genclone(y)



#expected vs observed heterozygosity might be interesting
# sumy <- summary(y)
# 
# png(filename = "heterozygosity.png", width=1000, height = 500)
# par(mfrow=c(1,2))
# hist(sumy$Hobs, main = "Histogram of observed heterozygosity")
# hist(sumy$Hexp, main = "Histogram of expected heterozygosity")
# dev.off()
# 
# hwp <- hw.test(y, B=0)
# 
# png(filename = "chi²-pval-hist.png", width=1000, height = 1000)
# hist(hwp[,3], breaks = 100, main = "chi²-test for HWE", xlab = "p-value")
# dev.off()


#snpStats: hard to import using GGtools
#ss <- vcf2sm('../data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf')



# plot HWE onto genome
# probably no interesting pattern
hw_df <- as.data.frame(hwp)
colnames(hw_df) <- c("chisq","df", "Pr_chisq")

vcf <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf")
hw_df$chr <- vcf@fix[,1]
hw_df$pos <- vcf@fix[,2]

p <- ggplot(hw_df, aes(pos,Pr_chisq)) + geom_point()

png(filename = "hwe_plotted.png", width=1000, height = 1000)
p + facet_wrap(vars(chr))
dev.off()

# some plots seem to be inverted compared to ggplot?
hwe1 <- subset(hw_df, chr=="chr01")
plot(hwe1$pos, hwe1$Pr_chisq)

hwe2 <- subset(hw_df, chr=="chr02")
plot(hwe2$pos, hwe2$Pr_chisq)

hwe11 <- subset(hw_df, chr=="chr11")
plot(hwe11$pos, hwe11$Pr_chisq)



# population strata needed?
# amova <- poppr.amova(y)

# only diploids
#filter_stats(y, plot = T)

# Index of Association
ia(y)

#long runtime
# res <- pair.ia(y)
# plot(res, low = "black", high = "green", index = "Ia")

# needs snplight
winia <- win.ia(y)
plot(winia, low = "black", high = "green", index = "Ia")


loctable <- locus_table(y)
#Hexp
hist(loctable[,3])

#Evenness?
hist(loctable[,4])

#1-D?
hist(loctable[,2])

# summary statistics
poppr(y)
