library(inbreedR)
library(vcfR)
library(reshape2)

vcf <- read.vcfR("../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_diploidized.vcf")


# extract genotypes
gt <- extract.gt(vcf)
# transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
# NA handling
gt[gt == "."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert
diploid_snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
check_data(diploid_snp_genotypes)

g2_dip_SNPs <- g2_snps(diploid_snp_genotypes, nperm = 100, nboot = 10, 
                       CI = 0.95, parallel = T, ncores = 11)

g2_dip_SNPs
plot(g2_dip_SNPs)
