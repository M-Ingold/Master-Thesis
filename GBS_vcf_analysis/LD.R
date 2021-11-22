library(VariantAnnotation)
library(updog)
library(ldsep)

chr2.gr <- GRanges(seqnames = "chr02", ranges=IRanges(start=1, end=100000000))
params <- ScanVcfParam(which=chr2.gr)
vcfGeno(object) <- "chr02"
vcf <- readVcf("data/VCF/freebayes_264_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf.gz", param = params)

geno(header(vcf))

sizemat <- geno(vcf)$DP
refmat <- geno(vcf)$RO
ploidy <- 4

mout <- multidog(refmat = refmat, 
                 sizemat = sizemat, 
                 ploidy = ploidy, 
                 model = "norm",
                 nc = 11)

plot(mout, indices = sample(1:nrow(vcf), 3))
msub <- filter_snp(x = mout, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95)
nrow(msub$snpdf)

varnames <- paste0("logL_", 0:ploidy)
varnames

larray <- format_multidog(x = msub, varname = varnames)
dim(larray)

pmmat <- format_multidog(x = msub, varname = "postmean")

like_ld <- mldest(geno = larray, K = ploidy, type = "comp", nc=11)
plot(like_ld)
par(mar = c(2.4, 2.8, 0, 0) + 0.5, mgp = c(1.8, 0.6, 0))
plot(mom_ld$r2, like_ld$r2, 
     xlab = expression(paste(textstyle(Naive), ~~hat(r)^2)), 
     ylab = expression(paste(textstyle(MLE), ~~hat(r)^2)), 
     pch  = 20)
abline(0, 1, lty = 2, col = 2)

ldmat <- format_lddf(obj = like_ld, element = "r2")
ldmat[1:4, 1:4]