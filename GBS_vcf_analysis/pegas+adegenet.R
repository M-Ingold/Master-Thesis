library(pegas)
#library(snpStats)
library(adegenet)
library(ggplot2)
#library(GGtools)

#x <- read.vcf('../data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf', to = 50000)

x <- read.vcf('../data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf', to = 7000)

X <- as.loci(x)

#convert to adegenet genind
y <- loci2genind(X, ploidy = 4)

#convert to genpop
z <- genind2genpop(y)

#expected vs observed heterozygosity might be interesting
sumy <- summary(y)

png(filename = "heterozygosity.png", width=1000, height = 500)
par(mfrow=c(1,2))
hist(sumy$Hobs)
hist(sumy$Hexp)
dev.off()

hwp <- hw.test(y, B=0)

png(filename = "chi²-pval-hist.png", width=1000, height = 1000)
hist(hwp[,3], breaks = 100, main = "chi²-test for HWE", xlab = "p-value")
dev.off()
#snpStats: hard to import using GGtools
#ss <- vcf2sm('../data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf')
